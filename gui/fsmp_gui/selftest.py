"""Health check of a Studio build: `FSMP-kMC Studio --selftest [report]`.

Verifies, without opening a window, that everything PyInstaller must carry
actually arrived: numpy, PySide6, RDKit with its MMFF parameter tables, the
bundled TMA model, and the engine binary next to the app. The release
workflow runs this on every platform right after assembling the bundle, so
a Studio that cannot type a molecule or find its engine never gets
published. Exit code 0 means healthy. The report goes to stdout and, when a
path is given, to that file: a windowed Windows build has no console, so
the workflow reads the file."""

import sys
import traceback


def run(report_path: str | None = None) -> int:
    lines = []
    ok = True
    try:
        import numpy
        lines.append(f"numpy {numpy.__version__}")
        import PySide6
        lines.append(f"PySide6 {PySide6.__version__}")
        import rdkit
        lines.append(f"rdkit {rdkit.__version__}")

        from .engine import app_root, find_engine
        from .generate import MMFFBackend
        from .mmff import mmff_pair_params
        from .molecule import Molecule

        model = app_root() / "samples" / "models" / "trimesic_acid.xyz"
        mol = Molecule.load_xyz(model)
        params = mmff_pair_params(mol)
        lines.append(f"MMFF typing of {model.name}: {len(params.elements)} "
                     f"atoms, charge sum {float(params.q.sum()):+.4f} e")

        # one real slab: proves RDKit's MMFF data tables were bundled and the
        # generator produces the validated hydrogen-bond well
        backend = MMFFBackend(mol)
        ang = numpy.deg2rad(numpy.arange(0, 361, 5))
        well = float(backend.slab(10.0, ang).min()) / 1000.0
        lines.append(f"pair well at r = 10 Å: {well:.1f} kJ/mol")
        if not -80.0 < well < -20.0:
            raise RuntimeError(f"implausible pair well: {well:.1f} kJ/mol")

        engine = find_engine()
        if engine is None:
            raise RuntimeError("engine binary not found next to the app")
        lines.append(f"engine: {engine[0]}")
    except Exception:
        ok = False
        lines.append(traceback.format_exc().rstrip())

    lines.append("SELFTEST OK" if ok else "SELFTEST FAILED")
    report = "\n".join(lines) + "\n"
    if report_path:
        with open(report_path, "w", encoding="utf-8") as f:
            f.write(report)
    try:
        sys.stdout.write(report)
    except Exception:
        pass                    # a windowed build may have no usable stdout
    return 0 if ok else 1
