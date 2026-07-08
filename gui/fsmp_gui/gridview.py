"""Shared 2D editor view: an angstrom grid with +y pointing up, wheel zoom
and a cursor-position signal. MoleculeCanvas and SiteCanvas build on it."""

import math

from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QPainter, QPen, QTransform
from PySide6.QtWidgets import QGraphicsScene, QGraphicsView

from . import theme


class GridView(QGraphicsView):
    cursorMoved = Signal(float, float)

    def __init__(self, parent=None):
        super().__init__(parent)
        scene = QGraphicsScene(self)
        # Editor items (bonds) are added and removed while an atom is dragged,
        # i.e. from inside the scene's own event delivery. The default BSP
        # index is not reentrant against that and crashes; with a handful of
        # items a linear scan (NoIndex) is free and safe.
        scene.setItemIndexMethod(QGraphicsScene.NoIndex)
        self.setScene(scene)
        self.scene().setSceneRect(-60, -60, 120, 120)
        self.setRenderHint(QPainter.Antialiasing)
        self.setDragMode(QGraphicsView.RubberBandDrag)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.setMouseTracking(True)
        self.setBackgroundBrush(QColor(theme.BG_PANEL))
        self.setFrameShape(QGraphicsView.NoFrame)
        self.setTransform(QTransform.fromScale(46, -46))

    def drawBackground(self, painter, rect):
        super().drawBackground(painter, rect)
        minor = QPen(QColor(theme.BORDER), 0)
        minor.setCosmetic(True)
        major = QPen(QColor(theme.BORDER).lighter(135), 0)
        major.setCosmetic(True)
        x0, x1 = math.floor(rect.left()), math.ceil(rect.right())
        y0, y1 = math.floor(rect.top()), math.ceil(rect.bottom())
        x = x0
        while x <= x1:
            painter.setPen(major if x % 5 == 0 else minor)
            painter.drawLine(QPointF(x, rect.top()), QPointF(x, rect.bottom()))
            x += 1
        y = y0
        while y <= y1:
            painter.setPen(major if y % 5 == 0 else minor)
            painter.drawLine(QPointF(rect.left(), y), QPointF(rect.right(), y))
            y += 1
        axes = QPen(QColor(theme.TEXT_DIM), 0)
        axes.setCosmetic(True)
        painter.setPen(axes)
        painter.drawLine(QPointF(rect.left(), 0), QPointF(rect.right(), 0))
        painter.drawLine(QPointF(0, rect.top()), QPointF(0, rect.bottom()))

    def mouseMoveEvent(self, event):
        pos = self.mapToScene(event.position().toPoint())
        self.cursorMoved.emit(pos.x(), pos.y())
        super().mouseMoveEvent(event)

    def wheelEvent(self, event):
        factor = 1.15 if event.angleDelta().y() > 0 else 1 / 1.15
        current = abs(self.transform().m11())
        if 4 < current * factor < 400:
            self.scale(factor, factor)

    def fit_points(self, points: list[tuple[float, float]], pad: float = 2.0):
        if points:
            xs = [p[0] for p in points]
            ys = [p[1] for p in points]
            rect = QRectF(QPointF(min(xs) - pad, min(ys) - pad),
                          QPointF(max(xs) + pad, max(ys) + pad))
        else:
            rect = QRectF(-6, -6, 12, 12)
        scale = min(self.viewport().width() / max(rect.width(), 1e-6),
                    self.viewport().height() / max(rect.height(), 1e-6))
        scale = min(max(scale, 4), 120)
        self.setTransform(QTransform.fromScale(scale, -scale))
        self.centerOn(rect.center())
