import argparse
import signal
import sys
import time
from pathlib import Path

import numpy as np

from fsmp_gui.generate import GridSpec, generate
from fsmp_gui.molecule import Atom, Molecule

from .backend import NNBackend
from .geometry import rotational_order, symmetry_error, symmetrize


def _say(text):
    print(text, file=sys.stderr, flush=True)


def _load(path):
    try:
        mol = Molecule.load_xyz(path)
    except (OSError, ValueError) as error:
        raise SystemExit(f"{path}: {error}")
    if not mol.atoms:
        raise SystemExit(f"{path}: the molecule has no atoms")
    return mol


def _coords(mol):
    return ([a.element for a in mol.atoms],
            np.array([(a.x, a.y) for a in mol.atoms], dtype=float))


def _rebuild(elements, xy, comment):
    return Molecule([Atom(e, float(xy[i, 0]), float(xy[i, 1]), 0.0)
                     for i, e in enumerate(elements)], comment)


def _calculator(args):
    from .aimnet import AIMNet2
    calc = AIMNet2(model=args.model, coulomb=args.coulomb, charge=args.charge)
    described = calc.describe()
    _say(f"model {described['model']} on {described['device']}, "
         f"short cutoff {described['short_cutoff']} A, "
         f"coulomb {described['coulomb']}, "
         f"declared cutoff {described['declared_cutoff']}")
    return calc


def _clock(seconds):
    if seconds < 90:
        return f"{seconds:.0f} s"
    if seconds < 5400:
        return f"{seconds / 60:.1f} min"
    return f"{seconds / 3600:.1f} h"


def cmd_optimize(args):
    from .relax import relax_planar

    mol = _load(args.molecule)
    elements, xy = _coords(mol)
    calc = _calculator(args)
    _say(f"relaxing {len(elements)} atoms, kept in the plane, fmax {args.fmax} eV/A")

    started = time.perf_counter()
    xy, converged, steps = relax_planar(elements, xy, calc, fmax=args.fmax,
                                        steps=args.max_steps)
    _say(f"{'converged' if converged else 'NOT converged'} after {steps} steps, "
         f"{_clock(time.perf_counter() - started)}")

    order = rotational_order(xy, elements, tol=args.symmetry_tol)
    error = symmetry_error(xy, elements, order)
    _say(f"rotational symmetry C{order}, off by {error:.4f} A")
    if args.symmetrize and order > 1:
        xy = symmetrize(xy, elements, order)
        _say(f"symmetrized to exact C{order}, off by {symmetry_error(xy, elements, order):.2e} A")

    out = Path(args.out)
    _rebuild(elements, xy, f"{mol.comment} relaxed with {args.model}").save_xyz(out)
    _say(f"wrote {out}")
    return 0


def cmd_grid(args):
    mol = _load(args.molecule)
    elements, xy = _coords(mol)

    order = rotational_order(xy, elements, tol=args.symmetry_tol)
    error = symmetry_error(xy, elements, order)
    if args.symmetrize and order > 1:
        xy = symmetrize(xy, elements, order)
        error = symmetry_error(xy, elements, order)
    _say(f"{len(elements)} atoms, rotational symmetry C{order}, off by {error:.2e} A")
    if order > 1 and error > args.symmetry_tol * 0.1:
        _say("note: folding averages the periods, which differ by that much; "
             "pass --symmetrize to make the molecule exactly symmetric first")

    fold = 360.0 / order if args.fold is None else args.fold
    radius = float(np.hypot(xy[:, 0], xy[:, 1]).max())
    r_min = round(radius, 2) if args.r_min is None else args.r_min

    calc = _calculator(args)
    charges = None
    if calc.cutoff is not None:
        from fsmp_gui.mmff import mmff_pair_params
        charges = np.asarray(mmff_pair_params(_rebuild(elements, xy, ""),
                                              int(args.charge)).q, dtype=float)
        _say(f"the model declares a cutoff, so the tail beyond {calc.cutoff} A "
             "comes from MMFF94 partial charges")

    backend = NNBackend(elements, xy, calc, args.da, charges=charges, order=order,
                        wall=args.wall, chunk=args.batch)
    spec = GridSpec(r_min=r_min, r_max=args.r_max, dr=args.dr, da=args.da,
                    fold_deg=fold, use_float=not args.double)

    per_row = len(backend._orbit(int(round(360.0 / order / args.da)))[0])
    total = spec.n_dist * per_row
    _say(f"grid r {r_min} to {args.r_max} A step {args.dr}, angles every {args.da} deg, "
         f"folded to {fold:.0f} deg")
    _say(f"{spec.n_dist} distances, {per_row} configurations each, {total} in total, "
         f"batch {backend.chunk}")

    stop = {"now": False}

    def on_sigint(signum, frame):
        stop["now"] = True
        _say("\nstopping after this distance")

    signal.signal(signal.SIGINT, on_sigint)
    started = time.perf_counter()

    def progress(done, count):
        elapsed = time.perf_counter() - started
        rate = done / elapsed if elapsed else 0.0
        left = (count - done) / rate if rate else 0.0
        print(f"\r  {done}/{count} distances, {_clock(elapsed)} elapsed, "
              f"{_clock(left)} left   ", end="", file=sys.stderr, flush=True)

    out = Path(args.out)
    done = generate(backend, spec, out, progress=progress, cancel=lambda: stop["now"])
    print(file=sys.stderr)
    elapsed = time.perf_counter() - started

    if not done:
        _say("cancelled, the partial file was removed")
        return 1

    evaluated = backend.stats["evaluated"]
    _say(f"wrote {out}, {out.stat().st_size / 1e6:.1f} MB in {_clock(elapsed)}")
    _say(f"{evaluated} configurations through the model "
         f"({evaluated / elapsed:.0f} per second), "
         f"{backend.stats['walled']} rejected by the {args.wall} A wall, "
         f"{backend.stats['rows_beyond_reach']} distances beyond its reach")
    return 0


def cmd_compare(args):
    from .compare import ComparisonError, compare, report
    try:
        result = compare(args.grid, args.reference, threshold=args.threshold)
    except ComparisonError as error:
        raise SystemExit(str(error))
    print(report(result, Path(args.grid).name, Path(args.reference).name))
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="python -m nnpot",
        description="Build an FSMP-kMC pair potential with a machine-learned model.")
    parser.add_argument("--model", default="aimnet2")
    parser.add_argument("--coulomb", default="simple", choices=("simple", "dsf", "ewald"))
    parser.add_argument("--charge", type=float, default=0.0)
    parser.add_argument("--symmetry-tol", type=float, default=0.05)
    sub = parser.add_subparsers(dest="command", required=True)

    opt = sub.add_parser("optimize", help="relax a molecule with the model, keeping it planar")
    opt.add_argument("--molecule", required=True)
    opt.add_argument("--out", required=True)
    opt.add_argument("--fmax", type=float, default=0.01)
    opt.add_argument("--max-steps", type=int, default=500)
    opt.add_argument("--symmetrize", action="store_true")
    opt.set_defaults(run=cmd_optimize)

    grid = sub.add_parser("grid", help="compute the pair potential and write it in the v2 format")
    grid.add_argument("--molecule", required=True)
    grid.add_argument("--out", required=True)
    grid.add_argument("--r-min", type=float, default=None)
    grid.add_argument("--r-max", type=float, default=30.0)
    grid.add_argument("--dr", type=float, default=0.02)
    grid.add_argument("--da", type=float, default=1.0)
    grid.add_argument("--fold", type=float, default=None)
    grid.add_argument("--wall", type=float, default=1.3)
    grid.add_argument("--batch", type=int, default=None)
    grid.add_argument("--double", action="store_true")
    grid.add_argument("--symmetrize", action="store_true")
    grid.set_defaults(run=cmd_grid)

    match = sub.add_parser("compare", help="compare a grid with a reference, point by point")
    match.add_argument("--grid", required=True)
    match.add_argument("--reference", required=True)
    match.add_argument("--threshold", type=float, default=5000.0)
    match.set_defaults(run=cmd_compare)

    args = parser.parse_args(argv)
    return args.run(args)


if __name__ == "__main__":
    sys.exit(main())
