"""__version__ must come from version.h, the single project-version source.

Run from the repository root:
    gui/.venv/Scripts/python gui/tests/test_version.py
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import fsmp_gui

REPO = Path(__file__).resolve().parents[2]


class TestVersion(unittest.TestCase):
    def test_parsed_from_version_h(self):
        self.assertNotEqual(fsmp_gui.__version__, "unknown")
        self.assertIn(f'"{fsmp_gui.__version__}"',
                      (REPO / "version.h").read_text(encoding="utf-8"))

    def test_looks_like_a_version(self):
        self.assertRegex(fsmp_gui.__version__, r"^\d+\.\d+\.\d+")


if __name__ == "__main__":
    unittest.main(verbosity=2)
