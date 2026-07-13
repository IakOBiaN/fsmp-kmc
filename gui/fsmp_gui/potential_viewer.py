"""A schematic view of two molecules for browsing a pair potential.

Molecule A sits at the origin turned by theta1, molecule B at distance r along
+x turned by theta2 - the same geometry the generator and the engine use. Each
molecule is drawn from a glyph (a list of coloured disks) with an orientation
arrow so the angle is visible even for symmetric shapes.
"""

import math

from PySide6.QtCore import Qt
from PySide6.QtGui import QBrush, QColor, QFont, QPainter, QPen, QTransform
from PySide6.QtWidgets import (QGraphicsEllipseItem, QGraphicsLineItem,
                               QGraphicsSimpleTextItem)

from . import theme
from .gridview import GridView

Glyph = list  # list[tuple[x, y, color_str, radius]]


class TwoMoleculeView(GridView):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._glyph = [(0.0, 0.0, theme.ACCENT, 0.5)]
        self._r, self._t1, self._t2 = 15.0, 0.0, 0.0

    def set_glyph(self, points) -> None:
        self._glyph = list(points) if points else [(0.0, 0.0, theme.ACCENT, 0.5)]
        self._redraw()

    def set_config(self, r: float, t1: float, t2: float) -> None:
        self._r, self._t1, self._t2 = r, t1, t2
        self._redraw()

    def _extent(self) -> float:
        return max((math.hypot(x, y) + rad for x, y, _, rad in self._glyph),
                   default=0.5)

    def base_distance(self) -> float:
        """A pleasant default separation, used when no potential is attached
        to provide the real distance range."""
        return round(max(2.4 * self._extent(), 6.0), 1)

    def _draw_molecule(self, cx: float, cy: float, deg: float, label: str) -> None:
        rad = math.radians(deg)
        c, s = math.cos(rad), math.sin(rad)
        for x, y, color, r in self._glyph:
            gx = cx + x * c - y * s
            gy = cy + x * s + y * c
            dot = QGraphicsEllipseItem(gx - r, gy - r, 2 * r, 2 * r)
            dot.setBrush(QBrush(QColor(color)))
            dot.setPen(QPen(QColor(color).darker(160), 0.04))
            dot.setZValue(1)
            self.scene().addItem(dot)
        reach = self._extent() + 0.8
        arrow = QGraphicsLineItem(cx, cy, cx + reach * c, cy + reach * s)
        arrow.setPen(QPen(QColor(theme.GOLD), 0.08, Qt.SolidLine, Qt.RoundCap))
        arrow.setZValue(2)
        self.scene().addItem(arrow)
        text = QGraphicsSimpleTextItem(label)
        font = QFont("Segoe UI", 10)
        font.setBold(True)
        text.setFont(font)
        text.setBrush(QBrush(QColor(theme.TEXT_DIM)))
        b = text.boundingRect()
        k = 0.9 / max(b.height(), 1e-6)
        text.setTransform(QTransform(k, 0, 0, -k, cx - k * b.width() / 2,
                                     cy - self._extent() - 1.0 + k * b.height() / 2))
        self.scene().addItem(text)

    def _redraw(self) -> None:
        self.scene().clear()
        axis = QGraphicsLineItem(0, 0, self._r, 0)
        pen = QPen(QColor(theme.TEXT_DIM), 0.05, Qt.DashLine)
        axis.setPen(pen)
        axis.setZValue(0)
        self.scene().addItem(axis)
        self._draw_molecule(0.0, 0.0, self._t1, "A")
        self._draw_molecule(self._r, 0.0, self._t2, "B")
        ext = self._extent()
        # the extra room below keeps the A/B labels inside the view
        self.fit_points([(-ext, -ext - 2.5), (-ext, ext),
                         (self._r + ext, ext),
                         (self._r + ext, -ext - 2.5)], pad=1.0)

    def showEvent(self, event):
        super().showEvent(event)
        self._redraw()
