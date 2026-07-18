"""Selection and rotation in the unit-cell editor: the canvas and the
placements table must point at the same molecule, the selected molecule
rotates by dragging its round handle, and the wheel always zooms (it must
never be overloaded: touchpad pinch arrives as Ctrl+wheel). Runs under the
offscreen Qt platform.

Run from the repository root:
    gui/.venv/Scripts/python gui/tests/test_unit_cell.py
"""

import math
import os
import sys
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from PySide6.QtCore import QEvent, QPoint, QPointF, Qt
from PySide6.QtGui import QMouseEvent, QWheelEvent
from PySide6.QtWidgets import QApplication

from fsmp_gui.placement_table import Placement, PlacementTableModel
from fsmp_gui.project import Project
from fsmp_gui.tabs.unit_cell_tab import UnitCellTab
from fsmp_gui.unit_cell_canvas import UnitCellCanvas

_app = QApplication.instance() or QApplication([])


def mouse(canvas, kind, scene_x, scene_y, modifiers=Qt.NoModifier):
    """A left-button mouse event aimed at a scene point."""
    vp = QPointF(canvas.mapFromScene(QPointF(scene_x, scene_y)))
    button = Qt.NoButton if kind == QEvent.Type.MouseMove else Qt.LeftButton
    return QMouseEvent(kind, vp, QPointF(canvas.mapToGlobal(vp.toPoint())),
                       button, Qt.LeftButton, modifiers)


def wheel(canvas, scene_x, scene_y, modifiers, delta=QPoint(0, 120)):
    """A wheel event aimed at a scene point, as the view would receive it."""
    vp = QPointF(canvas.mapFromScene(QPointF(scene_x, scene_y)))
    return QWheelEvent(vp, QPointF(canvas.mapToGlobal(vp.toPoint())), QPoint(),
                       delta, Qt.NoButton, modifiers, Qt.NoScrollPhase, False)


class TestSetPhi(unittest.TestCase):
    def test_normalizes_and_notifies(self):
        model = PlacementTableModel()
        model.set_placements([Placement(1.0, 2.0, 350.0)])
        got = []
        model.dataChanged.connect(
            lambda tl, br, roles=None: got.append((tl.row(), tl.column())))
        model.set_phi(0, 415.0)
        self.assertEqual(model.at(0).phi, 55.0)
        model.set_phi(0, -10.0)
        self.assertEqual(model.at(0).phi, 350.0)
        self.assertEqual(got, [(0, 2), (0, 2)])   # the phi column, twice


class TestHandleRotation(unittest.TestCase):
    def setUp(self):
        self.model = PlacementTableModel()
        self.canvas = UnitCellCanvas(self.model)
        self.canvas.resize(600, 400)
        self.canvas.set_cell(12.0, 12.0)
        self.model.set_placements([Placement(6.0, 6.0, 0.0),
                                   Placement(2.0, 2.0, 90.0)])
        self.picked = []
        self.canvas.moleculePicked.connect(self.picked.append)

    def grab_handle(self):
        hx, hy, _ = self.canvas._handle_geometry()
        self.canvas.mousePressEvent(
            mouse(self.canvas, QEvent.Type.MouseButtonPress, hx, hy))

    def test_handle_appears_only_with_a_selection(self):
        self.assertIsNone(self.canvas._handle_geometry())
        self.canvas.set_selected(0)
        hx, hy, knob = self.canvas._handle_geometry()
        # phi = 0: the handle sits along +x from the molecule centre
        self.assertGreater(hx, 6.0)
        self.assertAlmostEqual(hy, 6.0)
        self.assertGreater(knob, 0.0)

    def test_drag_rotates_the_selected_molecule(self):
        self.canvas.set_selected(0)
        self.grab_handle()
        self.assertEqual(self.canvas._rotate_row, 0)
        self.assertEqual(self.picked, [])       # grabbing keeps the selection
        self.canvas.mouseMoveEvent(
            mouse(self.canvas, QEvent.Type.MouseMove, 6.0, 9.0))
        self.assertAlmostEqual(self.model.at(0).phi, 90.0, delta=0.5)
        self.canvas.mouseReleaseEvent(
            mouse(self.canvas, QEvent.Type.MouseButtonRelease, 6.0, 9.0))
        self.assertIsNone(self.canvas._rotate_row)

    def test_shift_snaps_to_five_degrees(self):
        self.canvas.set_selected(0)
        self.grab_handle()
        rad = math.radians(47.0)
        self.canvas.mouseMoveEvent(
            mouse(self.canvas, QEvent.Type.MouseMove,
                  6.0 + 3.0 * math.cos(rad), 6.0 + 3.0 * math.sin(rad),
                  Qt.ShiftModifier))
        self.assertEqual(self.model.at(0).phi, 45.0)

    def test_click_on_a_molecule_picks_it(self):
        self.canvas.mousePressEvent(
            mouse(self.canvas, QEvent.Type.MouseButtonPress, 2.0, 2.0))
        self.assertEqual(self.picked, [1])
        self.assertEqual(self.canvas._drag_row, 1)

    def test_wheel_always_zooms(self):
        # Ctrl+wheel included: touchpads deliver pinch zoom exactly that way
        self.canvas.set_selected(0)
        for mods in (Qt.NoModifier, Qt.ControlModifier):
            before = abs(self.canvas.transform().m11())
            self.canvas.wheelEvent(wheel(self.canvas, 6.0, 6.0, mods))
            self.assertNotEqual(abs(self.canvas.transform().m11()), before,
                                msg=str(mods))
        self.assertEqual(self.model.at(0).phi, 0.0)   # nothing rotated


class TestSelectionSync(unittest.TestCase):
    """The tab wiring: one selection, visible on both sides."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(self._td.cleanup)
        project = Project.create(Path(self._td.name) / "proj", "t")
        self.tab = UnitCellTab(project)
        self.tab.model.set_placements([Placement(3.0, 3.0, 0.0),
                                       Placement(9.0, 9.0, 45.0)])

    def selected_rows(self):
        return [i.row()
                for i in self.tab.table.selectionModel().selectedRows()]

    def test_canvas_pick_selects_the_table_row(self):
        self.tab.canvas.moleculePicked.emit(1)
        self.assertEqual(self.selected_rows(), [1])
        self.assertEqual(self.tab.canvas._selected, 1)
        self.tab.canvas.moleculePicked.emit(-1)   # empty space clears
        self.assertEqual(self.selected_rows(), [])
        self.assertIsNone(self.tab.canvas._selected)

    def test_table_selection_highlights_the_canvas(self):
        self.tab.table.selectRow(0)
        self.assertEqual(self.tab.canvas._selected, 0)
        self.tab.table.clearSelection()
        self.assertIsNone(self.tab.canvas._selected)

    def test_removing_the_selected_row_clears_the_highlight(self):
        self.tab.table.selectRow(1)
        self.assertEqual(self.tab.canvas._selected, 1)
        self.tab.model.remove_rows([1])
        self.assertIsNone(self.tab.canvas._selected)


if __name__ == "__main__":
    unittest.main(verbosity=2)
