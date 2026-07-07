"""The site-model potential generator must reproduce the reference
forcefields/TMA_simple_2020.v2.bin from the bundled TMA site model.

Run from the repository root:
    gui/.venv/Scripts/python gui/tests/test_generate.py
"""

import array
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from fsmp_gui.generate import CAP_JMOL, GridSpec, _fold, _prepare, _slab
from fsmp_gui.sitemodel import SiteModel

REPO = Path(__file__).resolve().parents[2]
REFERENCE = REPO / "forcefields" / "TMA_simple_2020.v2.bin"
TMA = REPO / "models" / "TMA_simplified_2020.site"


class TestGenerator(unittest.TestCase):
    @unittest.skipUnless(REFERENCE.is_file(), "reference grid not present")
    def test_matches_reference_grid(self):
        model = SiteModel.load(TMA)
        lj, ch = _prepare(model)
        spec = GridSpec(7.6, 30.0, 0.02, 1.0, 120.0)
        ang = np.deg2rad(np.arange(int(round(360 / spec.da)) + 1) * spec.da)

        with REFERENCE.open("rb") as f:
            f.read(64)
            arr = array.array("f")
            arr.frombytes(f.read())
        n_ang = 121

        def ref(i, j, k):
            return arr[(i * n_ang + j) * n_ang + k]

        # a spread of distances and orientations, including near-core (capped),
        # attractive wells and the small-value cells where large terms cancel
        cells = [(1120, 0, 0), (1120, 30, 90), (500, 0, 60), (500, 60, 0),
                 (500, 45, 75), (160, 0, 60), (100, 0, 0), (100, 90, 15),
                 (5, 0, 60), (2, 20, 100)]
        cache = {}
        for i, j, k in cells:
            if i not in cache:
                r = 7.6 + i * 0.02
                slab = _slab(r, ang, lj, ch)
                np.clip(slab, -CAP_JMOL, CAP_JMOL, out=slab)
                cache[i] = _fold(slab, spec)
            got = cache[i][j, k]
            want = ref(i, j, k)
            # tolerance is the reference grid's own float32 / ASCII quantization
            self.assertLessEqual(abs(got - want), 0.1 + abs(want) * 3e-6,
                                 msg=f"cell {(i, j, k)}: {got} vs {want}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
