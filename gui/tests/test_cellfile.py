"""The .cell unit-cell format: round trip, error reporting, and the
committed reference cells extracted from StructureGenerator.h.

Run from the repository root:
    gui/.venv/Scripts/python gui/tests/test_cellfile.py
"""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsmp_gui.cellfile import load_cell, save_cell

REPO = Path(__file__).resolve().parents[2]
HCP = REPO / "cells" / "TMA_HCP_simple_2020.cell"
HEX = REPO / "cells" / "IPA_hexagones_qPBE_Dreiding_Dhb5.0.cell"


class TestRoundTrip(unittest.TestCase):
    def test_save_load(self):
        placements = [(0.0, 0.0, 90.0), (5.55, 9.6, 271.5)]
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "cell.cell"
            save_cell(path, 11.1, 19.2, placements, comment="test cell")
            cx, cy, back, comment = load_cell(path)
        self.assertEqual((cx, cy), (11.1, 19.2))
        self.assertEqual(comment, "test cell")
        for a, b in zip(back, placements):
            for u, v in zip(a, b):
                self.assertAlmostEqual(u, v, places=4)

    def test_errors(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "bad.cell"
            for text in ("", "x\n", "2\n11.1\n", "2\n11.1 19.2\n0 0 0\n",
                         "1\n11.1 19.2\n0 zero 0\n"):
                path.write_text(text)
                with self.assertRaises(ValueError, msg=repr(text)):
                    load_cell(path)


class TestReferenceCells(unittest.TestCase):
    @unittest.skipUnless(HCP.is_file(), "cells/ not present")
    def test_hcp_matches_the_engine_chain(self):
        cell_x, cell_y, placements, comment = load_cell(HCP)
        self.assertEqual((cell_x, cell_y), (11.1, 19.2))
        self.assertEqual(len(placements), 2)
        self.assertEqual(placements[0], (0.0, 0.0, 90.0))
        # 11.089 A at 59.967 deg from the first molecule
        self.assertAlmostEqual(placements[1][0], 5.5500, places=3)
        self.assertAlmostEqual(placements[1][1], 9.6002, places=3)
        self.assertEqual(placements[1][2], 90.0)
        self.assertIn("TMA_HCP_simple_2020", comment)

    @unittest.skipUnless(HEX.is_file(), "cells/ not present")
    def test_hexagones_count(self):
        _, _, placements, _ = load_cell(HEX)
        self.assertEqual(len(placements), 24)


if __name__ == "__main__":
    unittest.main(verbosity=2)
