"""Regression test for the molecule canvas bond bookkeeping.

Adding a third atom near two already-bonded atoms used to crash: rebuild()
calls scene.clear() (which deletes the bond line items), and decorate() then
called removeItem() on those now-deleted wrappers -> RuntimeError in the Qt
event loop. Runs under the offscreen Qt platform.

Run from the repository root:
    gui/.venv/Scripts/python gui/tests/test_canvas.py
"""

import os
import sys
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from PySide6.QtCore import QPointF
from PySide6.QtWidgets import QApplication

from fsmp_gui.atom_table import AtomTableModel
from fsmp_gui.canvas import MoleculeCanvas
from fsmp_gui.molecule import Atom

_app = QApplication.instance() or QApplication([])


class TestMoleculeCanvas(unittest.TestCase):
    def setUp(self):
        self.model = AtomTableModel()
        self.canvas = MoleculeCanvas(self.model)

    def test_incremental_add_with_bonds(self):
        # each add triggers rowsInserted -> rebuild -> scene.clear()
        for el, x, y in [("C", 0, 0), ("O", 1.2, 0), ("N", 0, 1.2),
                         ("S", 1.2, 1.2), ("F", -1.1, 0)]:
            self.model.add_atom(Atom(el, x, y))
        self.assertEqual(len(self.model.molecule.atoms), 5)
        # bonds were rebuilt without touching stale line items
        self.assertGreater(len(self.canvas._bonds), 0)

    def test_remove_all(self):
        for el, x, y in [("C", 0, 0), ("O", 1.2, 0), ("N", 0, 1.2)]:
            self.model.add_atom(Atom(el, x, y))
        while self.model.molecule.atoms:
            self.model.remove_rows([0])
        self.assertEqual(self.canvas._bonds, [])

    def test_drag_repositions_bonds(self):
        for el, x, y in [("C", 0, 0), ("O", 1.2, 0), ("N", 0, 1.2)]:
            self.model.add_atom(Atom(el, x, y))
        before = len(self.canvas._bonds)
        # simulate a drag of atom 0: setPos fires itemChange -> item_dragged
        for i in range(50):
            self.canvas._items[0].setPos(QPointF(0.02 * i, 0.0))
        # dragging must not leave stale/duplicated lines
        self.assertLessEqual(len(self.canvas._bonds), before + 2)
        self.assertAlmostEqual(self.model.molecule.atoms[0].x, 0.98, places=2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
