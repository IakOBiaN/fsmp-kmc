"""Engine-integration helpers: the chained-polar unit_cell parameter, the
xyz animation decoding and the final report parsing must agree with the
engine's conventions. The last test runs the real engine on the committed
test grid (skipped when the built engine or the grid is absent).

Run from the repository root:
    gui/.venv/Scripts/python gui/tests/test_engine.py
"""

import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsmp_gui.engine import (ANIMATION, PARAMETERS, FrameReader, app_root,
                             chain_to_placements, decode_frame, final_energy,
                             final_params, find_engine, parameter_text,
                             placements_to_chain)

REPO = Path(__file__).resolve().parents[2]
GRID = REPO / "tests" / "data" / "TMA_simple_2020_s4.v2.bin"
MODEL = REPO / "samples" / "models" / "trimesic_acid.xyz"


def rotated(offsets, x, y, phi):
    rad = math.radians(phi)
    c, s = math.cos(rad), math.sin(rad)
    return [(x + ox * c - oy * s, y + ox * s + oy * c) for ox, oy in offsets]


def frame_text(offsets, cell_x, cell_y, placements):
    """A frame as the engine writes it: the 3x3 tiling of the cell, the
    first-tile molecules first."""
    atoms = []
    for ix in range(3):
        for iy in range(3):
            for x, y, phi in placements:
                for px, py in rotated(offsets, x + ix * cell_x,
                                      y + iy * cell_y, phi):
                    atoms.append(f"C {px:.8f} {py:.8f} 0 1 0")
    head = [str(len(atoms)),
            f'Lattice="{cell_x * 3} 0 0 0 {cell_y * 3} 0 0 0 1" '
            'Properties=species:S:1:pos:R:3:lambda:R:1:u_ext:R:1']
    return "\n".join(head + atoms) + "\n"


class Helpers(unittest.TestCase):
    def assert_circular(self, a, b, period, tol, msg=""):
        d = abs(a - b) % period
        self.assertLessEqual(min(d, period - d), tol,
                             msg=f"{msg}: {a} vs {b} (period {period})")

    def assert_placements(self, got, want, cell_x, cell_y, tol=1e-6,
                          phi_tol=1e-6):
        self.assertEqual(len(got), len(want))
        for g, w in zip(got, want):
            self.assert_circular(g[0], w[0], cell_x, tol, "x")
            self.assert_circular(g[1], w[1], cell_y, tol, "y")
            self.assert_circular(g[2], w[2], 360.0, phi_tol, "phi")


class TestChain(Helpers):
    def test_round_trip(self):
        cell_x, cell_y = 12.5, 15.0
        placements = [(1.2, 3.4, 10.0), (7.0, 2.0, 250.0), (0.5, 14.0, 84.0)]
        chain = placements_to_chain(cell_x, cell_y, placements)
        values = [float(t) for t in chain.split()]
        self.assertEqual(len(values), 3 + 3 * len(placements))
        cx, cy, back = chain_to_placements(values)
        self.assertAlmostEqual(cx, cell_x)
        self.assertAlmostEqual(cy, cell_y)
        self.assert_placements(back, placements, cell_x, cell_y, tol=1e-5,
                               phi_tol=1e-5)

    def test_first_molecule_at_origin(self):
        chain = placements_to_chain(10.0, 10.0, [(0.0, 0.0, 84.0)])
        values = [float(t) for t in chain.split()]
        _, _, back = chain_to_placements(values)
        self.assert_placements(back, [(0.0, 0.0, 84.0)], 10.0, 10.0)


class TestFrameDecode(Helpers):
    OFFSETS = [(0.0, 0.0), (1.5, 0.0), (0.0, 2.0)]

    def test_two_molecules(self):
        cell_x, cell_y = 12.0, 16.0
        placements = [(3.0, 4.0, 30.0), (8.0, 10.0, 300.0)]
        lines = frame_text(self.OFFSETS, cell_x, cell_y,
                           placements).splitlines()
        cx, cy, got = decode_frame(lines, self.OFFSETS, 2)
        self.assertAlmostEqual(cx, cell_x)
        self.assertAlmostEqual(cy, cell_y)
        self.assert_placements(got, placements, cell_x, cell_y, tol=1e-6,
                               phi_tol=1e-6)

    def test_single_site_has_no_orientation(self):
        lines = frame_text([(0.0, 0.0)], 10.0, 10.0,
                           [(2.0, 3.0, 45.0)]).splitlines()
        _, _, got = decode_frame(lines, [(0.0, 0.0)], 1)
        self.assert_placements(got, [(2.0, 3.0, 0.0)], 10.0, 10.0)


class TestFrameReader(Helpers):
    def test_incremental_frames(self):
        offsets = [(0.0, 0.0), (1.0, 0.0)]
        first = frame_text(offsets, 10.0, 10.0, [(2.0, 2.0, 0.0)])
        second = frame_text(offsets, 10.0, 10.0, [(4.0, 5.0, 90.0)])
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / ANIMATION
            reader = FrameReader(path, offsets, 1)
            self.assertIsNone(reader.poll())          # no file yet
            cut = len(first) // 2
            path.write_bytes(first[:cut].encode())
            self.assertIsNone(reader.poll())          # half a frame
            with open(path, "ab") as f:
                f.write((first[cut:] + second).encode())
            frame = reader.poll()                     # the newest wins
            self.assertIsNotNone(frame)
            _, _, got = frame
            self.assert_placements(got, [(4.0, 5.0, 90.0)], 10.0, 10.0)
            self.assertIsNone(reader.poll())          # nothing new


class TestFinalReport(unittest.TestCase):
    OUTPUT = ("Cell scaling: factor 1.05, energy -50 kJ/mol per molecule\n"
              "Density: 1.5\t Energy: -60.1\n"
              "Optimization converged after 100 iterations (10 accepted)\n"
              "Final energy per molecule: -62.2276 kJ/mol\n"
              "Final density: 1.6656 mkmol/m2\n"
              "Final params: \n"
              "Number 0: 2\nNumber 1: 11.1\nNumber 2: 19.2\n"
              "Number 3: 0\nNumber 4: 0\nNumber 5: 84.5\n"
              "Number 6: 11.089\nNumber 7: 59.967\nNumber 8: 96.2\n")

    def test_parse(self):
        params = final_params(self.OUTPUT)
        self.assertEqual(params, [2, 11.1, 19.2, 0, 0, 84.5,
                                  11.089, 59.967, 96.2])
        self.assertAlmostEqual(final_energy(self.OUTPUT), -62.2276)

    def test_absent(self):
        self.assertIsNone(final_params("no report here"))
        self.assertIsNone(final_energy("no report here"))


class TestFindEngine(unittest.TestCase):
    """The resolution order. The frozen-app branch is faked by setting
    sys.frozen the way PyInstaller does."""

    def setUp(self):
        self._saved = os.environ.pop("FSMP_ENGINE", None)

    def tearDown(self):
        if self._saved is not None:
            os.environ["FSMP_ENGINE"] = self._saved

    def test_env_override_wins(self):
        with tempfile.TemporaryDirectory() as td:
            fake = Path(td) / "custom-engine.exe"
            fake.write_bytes(b"")
            os.environ["FSMP_ENGINE"] = str(fake)
            self.assertEqual(find_engine(), [str(fake)])

    def test_missing_override_is_ignored(self):
        os.environ["FSMP_ENGINE"] = str(Path("nowhere") / "fsmp.exe")
        command = find_engine()
        if command is not None:
            self.assertNotIn("nowhere", command[0])

    def test_frozen_app_looks_next_to_the_executable(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td).resolve()
            name = "fsmp.exe" if os.name == "nt" else "fsmp"
            engine = root / name
            engine.write_bytes(b"")
            old_exe = sys.executable
            sys.frozen = True
            sys.executable = str(root / "FSMP-kMC Studio.exe")
            try:
                self.assertEqual(app_root(), root)
                self.assertEqual(find_engine(), [str(engine)])
            finally:
                del sys.frozen
                sys.executable = old_exe

    def test_mac_bundle_data_sits_next_to_the_app(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td).resolve()
            macos = root / "FSMP-kMC Studio.app" / "Contents" / "MacOS"
            macos.mkdir(parents=True)
            old_exe = sys.executable
            sys.frozen = True
            sys.executable = str(macos / "FSMP-kMC Studio")
            try:
                self.assertEqual(app_root(), root)
            finally:
                del sys.frozen
                sys.executable = old_exe


@unittest.skipUnless(find_engine() and GRID.is_file() and MODEL.is_file(),
                     "engine or test grid not present")
class TestEndToEnd(Helpers):
    def test_optimize_one_molecule(self):
        lines = MODEL.read_text(encoding="utf-8").splitlines()
        n = int(lines[0].split()[0])
        rows = [ln.split() for ln in lines[2:2 + n]]
        offsets = [(float(r[1]), float(r[2])) for r in rows]
        # ignore_cleanup_errors: on Windows the directory handle can stay
        # held a moment after the engine exits
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
            run = Path(td)
            model = [str(n), "test"] + [f"{r[0]} {r[1]} {r[2]} 0.0"
                                        for r in rows]
            (run / "model.xyz").write_text("\n".join(model) + "\n")
            rel = os.path.relpath(GRID, run).replace("\\", "/")
            chain = placements_to_chain(11.6, 19.8, [(0.0, 0.0, 84.0)])
            (run / PARAMETERS).write_text(parameter_text(rel, chain))
            proc = subprocess.run(find_engine() + [PARAMETERS], cwd=run,
                                  capture_output=True, text=True, timeout=600)
            self.assertEqual(proc.returncode, 0,
                             msg=proc.stdout[-2000:] + proc.stderr[-2000:])

            params = final_params(proc.stdout)
            self.assertIsNotNone(params, msg=proc.stdout[-2000:])
            cell_x, cell_y, final = chain_to_placements(params)
            energy = final_energy(proc.stdout)
            self.assertIsNotNone(energy)
            self.assertLess(energy, 0.0)

            # the last animation frame must decode to the reported cell
            reader = FrameReader(run / ANIMATION, offsets, 1)
            frame = reader.poll()
            self.assertIsNotNone(frame)
            fx, fy, got = frame
            self.assertAlmostEqual(fx, cell_x, delta=1e-3)
            self.assertAlmostEqual(fy, cell_y, delta=1e-3)
            self.assert_placements(got, final, cell_x, cell_y, tol=0.02,
                                   phi_tol=0.2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
