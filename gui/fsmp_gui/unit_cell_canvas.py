"""Editor canvas for the rough unit cell: a rectangular cell with molecule
copies that can be placed, dragged and (via the table) rotated. Faint periodic
images show how the copies pack across cell boundaries.

Dragging is handled manually (hit-test + full redraw) rather than via movable
Qt items: a molecule is many disks, and the scene uses NoIndex so redrawing on
every mouse move is cheap and safe."""

import math

from PySide6.QtCore import Qt
from PySide6.QtGui import QBrush, QColor, QPainter, QPen
from PySide6.QtWidgets import QGraphicsEllipseItem, QGraphicsRectItem

from . import theme
from .canvas import Mode
from .glyph import fallback_glyph
from .gridview import GridView
from .placement_table import Placement, PlacementTableModel


class UnitCellCanvas(GridView):
    def __init__(self, model: PlacementTableModel, parent=None):
        super().__init__(parent)
        self.model = model
        self.glyph = fallback_glyph()
        self.cell_x, self.cell_y = 10.0, 10.0
        self.mode = Mode.SELECT
        self._drag_row: int | None = None

        model.modelReset.connect(self._redraw_fit)
        model.rowsInserted.connect(self._redraw_fit)
        model.rowsRemoved.connect(self._redraw_fit)
        model.dataChanged.connect(self._redraw)
        self._redraw()

    # -- configuration -----------------------------------------------------

    def set_glyph(self, points) -> None:
        self.glyph = list(points) if points else fallback_glyph()
        self._redraw_fit()

    def set_cell(self, cell_x: float, cell_y: float) -> None:
        self.cell_x, self.cell_y = cell_x, cell_y
        self._redraw_fit()

    def set_mode(self, mode: Mode) -> None:
        self.mode = mode

    def _glyph_extent(self) -> float:
        return max((math.hypot(x, y) + r for x, y, _, r in self.glyph), default=0.5)

    # -- drawing -----------------------------------------------------------

    def _draw_molecule(self, px: float, py: float, phi: float, opacity: float):
        rad = math.radians(phi)
        c, s = math.cos(rad), math.sin(rad)
        for x, y, color, r in self.glyph:
            gx = px + x * c - y * s
            gy = py + x * s + y * c
            dot = QGraphicsEllipseItem(gx - r, gy - r, 2 * r, 2 * r)
            col = QColor(color)
            dot.setBrush(QBrush(col))
            dot.setPen(QPen(col.darker(160), 0.04))
            dot.setOpacity(opacity)
            dot.setZValue(1 if opacity >= 1.0 else 0)
            self.scene().addItem(dot)

    def _redraw(self) -> None:
        self.scene().clear()
        rect = QGraphicsRectItem(0, 0, self.cell_x, self.cell_y)
        rect.setPen(QPen(QColor(theme.ACCENT), 0.06))
        rect.setBrush(QBrush(QColor(theme.ACCENT_BG)))
        rect.setZValue(-1)
        self.scene().addItem(rect)
        shifts = (-1, 0, 1)
        for p in self.model.placements:
            for ix in shifts:
                for iy in shifts:
                    if ix == 0 and iy == 0:
                        continue
                    self._draw_molecule(p.x + ix * self.cell_x,
                                        p.y + iy * self.cell_y, p.phi, 0.22)
        for p in self.model.placements:
            self._draw_molecule(p.x, p.y, p.phi, 1.0)

    def _redraw_fit(self, *args) -> None:
        self._redraw()
        m = self._glyph_extent()
        self.fit_points([(-m, -m), (self.cell_x + m, self.cell_y + m)], pad=1.5)

    # -- manual interaction ------------------------------------------------

    def _nearest_row(self, sx: float, sy: float) -> int | None:
        best, best_d = None, self._glyph_extent() + 0.5
        for i, p in enumerate(self.model.placements):
            d = math.hypot(p.x - sx, p.y - sy)
            if d < best_d:
                best, best_d = i, d
        return best

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            pos = self.mapToScene(event.position().toPoint())
            if self.mode == Mode.ADD:
                self.model.add(Placement(round(pos.x(), 3), round(pos.y(), 3), 0.0))
                return
            row = self._nearest_row(pos.x(), pos.y())
            if self.mode == Mode.DELETE and row is not None:
                self.model.remove_rows([row])
                return
            if self.mode == Mode.SELECT and row is not None:
                self._drag_row = row
                return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        pos = self.mapToScene(event.position().toPoint())
        self.cursorMoved.emit(pos.x(), pos.y())
        if self._drag_row is not None:
            self.model.move(self._drag_row, round(pos.x(), 3), round(pos.y(), 3))
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        self._drag_row = None
        super().mouseReleaseEvent(event)

    def paintEvent(self, event):
        # antialiased dots
        self.setRenderHint(QPainter.Antialiasing, True)
        super().paintEvent(event)
