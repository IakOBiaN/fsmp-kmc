import argparse
import shutil
import signal
import struct
import sys
import time
from pathlib import Path

import numpy as np

from fsmp_gui.generate import HEADER_BYTES, MAGIC, VERSION, GridSpec, generate
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


def _absorb_chunk(partial, chunk, per_angle, row_bytes, spec, args):
    if not chunk.exists():
        return
    if not partial.exists():
        raise SystemExit(f"{chunk} continues a file that is gone; delete it to start over")
    raw = chunk.read_bytes()[:HEADER_BYTES]
    if len(raw) < HEADER_BYTES or raw[:4] != MAGIC:
        raise SystemExit(f"{chunk} is not a v2 grid; delete it to start over")
    version, dtype, n_dist, n_ang = struct.unpack("<4I", raw[4:20])
    min_dist, dr, da, fold = struct.unpack("<4d", raw[20:52])
    done = (partial.stat().st_size - HEADER_BYTES) // row_bytes
    follows = (version == VERSION and n_ang == per_angle
               and dtype == (0 if args.double else 1)
               and abs(min_dist - (spec.r_min + done * spec.dr)) < 1e-9
               and abs(dr - spec.dr) < 1e-12
               and abs(da - spec.da) < 1e-12
               and abs(fold - spec.fold_deg) < 1e-9)
    if not follows:
        raise SystemExit(f"{chunk} does not continue {partial.name}; "
                         "delete it to start over")
    rows = (chunk.stat().st_size - HEADER_BYTES) // row_bytes
    if rows:
        with partial.open("ab") as whole, chunk.open("rb") as rest:
            rest.seek(HEADER_BYTES)
            whole.write(rest.read(rows * row_bytes))
    chunk.unlink()


def _resume_point(partial, spec, per_angle, row_bytes, args):
    if not partial.exists():
        return 0
    raw = partial.read_bytes()[:HEADER_BYTES]
    if len(raw) < HEADER_BYTES or raw[:4] != MAGIC:
        raise SystemExit(f"{partial} is not a v2 grid; delete it to start over")
    version, dtype, n_dist, n_ang = struct.unpack("<4I", raw[4:20])
    min_dist, dr, da, fold = struct.unpack("<4d", raw[20:52])
    same = (version == VERSION and n_ang == per_angle and n_dist == spec.n_dist
            and dtype == (0 if args.double else 1)
            and abs(min_dist - spec.r_min) < 1e-9
            and abs(dr - spec.dr) < 1e-12
            and abs(da - spec.da) < 1e-12
            and abs(fold - spec.fold_deg) < 1e-9)
    if not same:
        raise SystemExit(f"{partial} was written with different settings; "
                         "delete it to start over")
    ready = (partial.stat().st_size - HEADER_BYTES) // row_bytes
    keep = HEADER_BYTES + ready * row_bytes
    if partial.stat().st_size != keep:
        with partial.open("r+b") as f:
            f.truncate(keep)
    return int(min(ready, spec.n_dist))


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

    per_angle = int(round(fold / args.da)) + 1
    per_row = len(backend._orbit(int(round(360.0 / order / args.da)))[0])
    total = spec.n_dist * per_row
    _say(f"grid r {r_min} to {args.r_max} A step {args.dr}, angles every {args.da} deg, "
         f"folded to {fold:.0f} deg")
    _say(f"{spec.n_dist} distances, {per_row} configurations each, {total} in total, "
         f"batch {backend.chunk}")

    out = Path(args.out)
    partial = out.with_name(out.name + ".partial")
    row_bytes = per_angle * per_angle * (8 if args.double else 4)
    chunk = partial.with_name(partial.name + ".chunk")
    _absorb_chunk(partial, chunk, per_angle, row_bytes, spec, args)
    ready = _resume_point(partial, spec, per_angle, row_bytes, args)
    if ready:
        _say(f"resuming: {ready} of {spec.n_dist} distances are already on disk")
        spec = GridSpec(r_min=r_min + ready * args.dr, r_max=args.r_max, dr=args.dr,
                        da=args.da, fold_deg=fold, use_float=not args.double)
        target = chunk
    else:
        target = partial

    stop = {"now": False}

    def on_sigint(signum, frame):
        stop["now"] = True
        _say("\nstopping after this distance")

    signal.signal(signal.SIGINT, on_sigint)
    started = time.perf_counter()

    def progress(count, of):
        elapsed = time.perf_counter() - started
        rate = count / elapsed if elapsed else 0.0
        left = (of - count) / rate if rate else 0.0
        print(f"\r  {ready + count}/{ready + of} distances, {_clock(elapsed)} elapsed, "
              f"{_clock(left)} left   ", end="", file=sys.stderr, flush=True)

    finished = generate(backend, spec, target, progress=progress,
                        cancel=lambda: stop["now"])
    print(file=sys.stderr)
    elapsed = time.perf_counter() - started

    if not finished:
        _say(f"stopped; {partial} keeps the finished distances, "
             "run the same command again to continue")
        return 1

    if target != partial:
        with partial.open("ab") as whole, target.open("rb") as rest:
            rest.seek(HEADER_BYTES)
            shutil.copyfileobj(rest, whole)
        target.unlink()
    partial.replace(out)

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
