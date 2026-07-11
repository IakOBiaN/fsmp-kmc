"""The simulation-cell preview must mirror the engine exactly:
generate_elongated_cell (StructureGenerator.h) for the geometry and
fields.h for the damping and external-field profiles.

Run from the repository root:
    gui/.venv/Scripts/python gui/tests/test_simcell.py
"""

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsmp_gui.simcell import (SimCellSpec, damping_sqrt, external_u,
                              lambda_profile, tile_unit_cell, u_profile)

UC = {"cell_x": 11.1, "cell_y": 19.2,
      "molecules": [{"x": 0.0, "y": 0.0, "phi": 90.0},
                    {"x": 5.55, "y": 9.6, "phi": 90.0}]}


class TestFields(unittest.TestCase):
    def test_damping_pins(self):
        l0, lm = 0.5, 0.1
        self.assertEqual(damping_sqrt(0.0, l0, lm), 1.0)
        # ksi = 6: l0 + (1-l0)*(7-6)^2*(6-4)/4 = l0 + (1-l0)/2
        self.assertAlmostEqual(damping_sqrt(6.0, l0, lm), 0.75)
        self.assertEqual(damping_sqrt(8.0, l0, lm), l0)
        # ksi = 10: lm + (l0-lm)*(11-10)^2*(10-8)/4 = lm + (l0-lm)/2
        self.assertAlmostEqual(damping_sqrt(10.0, l0, lm), 0.3)
        self.assertEqual(damping_sqrt(12.0, l0, lm), lm)

    def test_damping_is_continuous(self):
        l0, lm = 0.577, 0.2
        for edge in (5.0, 7.0, 9.0, 11.0):
            a = damping_sqrt(edge - 1e-9, l0, lm)
            b = damping_sqrt(edge + 1e-9, l0, lm)
            self.assertAlmostEqual(a, b, places=6, msg=f"jump at ksi={edge}")

    def test_external_pins(self):
        um = -5000.0
        self.assertEqual(external_u(0.0, um), 0.0)
        self.assertAlmostEqual(external_u(6.0, um), um / 2.0)
        self.assertEqual(external_u(8.0, um), um)
        for edge in (5.0, 7.0):
            self.assertAlmostEqual(external_u(edge - 1e-9, um),
                                   external_u(edge + 1e-9, um), places=4)

    def test_profiles_in_x(self):
        spec = SimCellSpec(temp_transition=900.0)
        lx = 400.0
        # center: full interactions; plateau: T/T_tz; edge: lambdam^2
        self.assertEqual(lambda_profile(lx / 2, lx, 300.0, spec), 1.0)
        x_plateau = lx / 2 + 8.0 * lx / 32.0
        self.assertAlmostEqual(lambda_profile(x_plateau, lx, 300.0, spec),
                               300.0 / 900.0)
        self.assertAlmostEqual(lambda_profile(0.0, lx, 300.0, spec),
                               spec.lambdam ** 2)
        self.assertEqual(u_profile(lx / 2, lx, -5000.0), 0.0)
        self.assertEqual(u_profile(0.0, lx, -5000.0), -5000.0)


class TestLayout(unittest.TestCase):
    def test_geometry_matches_engine_formula(self):
        spec = SimCellSpec(uc_in_x=22, uc_in_y=6, free_space=0.24)
        layout, placements = tile_unit_cell(UC, spec)
        self.assertAlmostEqual(layout.slab_x, 11.1 * 22)
        self.assertAlmostEqual(layout.lx, 11.1 * 22 / (1 - 2 * 0.24))
        self.assertAlmostEqual(layout.ly, 19.2 * 6)
        self.assertEqual(layout.n_molecules, 2 * 22 * 6)
        self.assertEqual(len(placements), layout.n_molecules)
        # the slab is centered and every molecule stays inside it
        self.assertAlmostEqual(layout.slab_x0,
                               (layout.lx - layout.slab_x) / 2.0)
        xs = [p[0] for p in placements]
        self.assertGreaterEqual(min(xs), layout.slab_x0)
        self.assertLessEqual(max(xs), layout.slab_x0 + layout.slab_x)

    def test_zone_positions_are_symmetric(self):
        layout, _ = tile_unit_cell(UC, SimCellSpec())
        for ksi in (5.0, 7.0, 9.0, 11.0):
            left, right = layout.zone_x(ksi)
            self.assertAlmostEqual(layout.lx - right, left)
            self.assertAlmostEqual(32.0 * abs(right - layout.lx / 2) / layout.lx,
                                   ksi)


if __name__ == "__main__":
    unittest.main(verbosity=2)
