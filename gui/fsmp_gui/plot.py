"""A small dependency-free XY plot: line + markers, autoscale, nice ticks.
Enough for the run statistics (any column against any column) and cheap to
redraw live while the statistics file grows."""

import math

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QFontMetricsF, QPainter, QPen
from PySide6.QtWidgets import QWidget

from . import theme


def nice_ticks(lo: float, hi: float, target: int = 5) -> list:
    """Tick positions covering [lo, hi] with a 1/2/5 x 10^k step."""
    if not (math.isfinite(lo) and math.isfinite(hi)) or hi <= lo:
        return [lo]
    raw = (hi - lo) / max(target, 2)
    power = 10.0 ** math.floor(math.log10(raw))
    for mult in (1.0, 2.0, 5.0, 10.0):
        step = mult * power
        if raw <= step:
            break
    first = math.ceil(lo / step) * step
    ticks = []
    t = first
    while t <= hi + step * 1e-9:
        ticks.append(round(t, 12))
        t += step
    return ticks


class PlotWidget(QWidget):
    # the top and bottom margins leave a full text line for the axis
    # titles, so they never collide with the tick labels
    MARGIN_L, MARGIN_R, MARGIN_T, MARGIN_B = 74, 18, 34, 60

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(220)
        self._x, self._y = [], []
        self._xlabel, self._ylabel = "", ""

    def set_data(self, x: list, y: list, xlabel: str, ylabel: str) -> None:
        # keep only finite pairs; the engine writes nan for empty gas phases
        pairs = [(a, b) for a, b in zip(x, y)
                 if math.isfinite(a) and math.isfinite(b)]
        self._x = [p[0] for p in pairs]
        self._y = [p[1] for p in pairs]
        self._xlabel, self._ylabel = xlabel, ylabel
        self.update()

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        p.fillRect(self.rect(), QColor(theme.BG_PANEL))
        rect = self.rect().adjusted(self.MARGIN_L, self.MARGIN_T,
                                    -self.MARGIN_R, -self.MARGIN_B)
        if not self._x:
            p.setPen(QColor(theme.TEXT_DIM))
            p.drawText(self.rect(), Qt.AlignCenter, "no data points yet")
            return

        x_lo, x_hi = min(self._x), max(self._x)
        y_lo, y_hi = min(self._y), max(self._y)
        x_pad = (x_hi - x_lo) * 0.06 or max(abs(x_lo), 1.0) * 0.05
        y_pad = (y_hi - y_lo) * 0.08 or max(abs(y_lo), 1.0) * 0.05
        x_lo, x_hi = x_lo - x_pad, x_hi + x_pad
        y_lo, y_hi = y_lo - y_pad, y_hi + y_pad
        sx = lambda v: rect.left() + (v - x_lo) / (x_hi - x_lo) * rect.width()
        sy = lambda v: rect.bottom() - (v - y_lo) / (y_hi - y_lo) * rect.height()

        grid = QPen(QColor(theme.BORDER), 1)
        text = QColor(theme.TEXT_DIM)
        fm = QFontMetricsF(p.font())
        for t in nice_ticks(x_lo, x_hi):
            px = sx(t)
            p.setPen(grid)
            p.drawLine(QPointF(px, rect.top()), QPointF(px, rect.bottom()))
            p.setPen(text)
            label = f"{t:g}"
            p.drawText(QPointF(px - fm.horizontalAdvance(label) / 2,
                               rect.bottom() + fm.height() + 4), label)
        for t in nice_ticks(y_lo, y_hi):
            py = sy(t)
            p.setPen(grid)
            p.drawLine(QPointF(rect.left(), py), QPointF(rect.right(), py))
            p.setPen(text)
            label = f"{t:g}"
            p.drawText(QPointF(rect.left() - fm.horizontalAdvance(label) - 8,
                               py + fm.height() / 3), label)
        p.setPen(QPen(QColor(theme.BORDER).lighter(135), 1))
        p.drawRect(rect)

        # axis titles: x centered below, y along the top-left corner
        p.setPen(text)
        p.drawText(QPointF(rect.center().x()
                           - fm.horizontalAdvance(self._xlabel) / 2,
                           self.height() - 8), self._xlabel)
        p.drawText(QPointF(6, 6 + fm.ascent()), self._ylabel)

        # the series: points in file order, connected when there are several
        points = [QPointF(sx(a), sy(b)) for a, b in zip(self._x, self._y)]
        if len(points) > 1:
            pen = QPen(QColor(theme.ACCENT_HOVER), 1.6)
            p.setPen(pen)
            for a, b in zip(points, points[1:]):
                p.drawLine(a, b)
        p.setPen(QPen(QColor(theme.GOLD), 1))
        p.setBrush(QColor(theme.GOLD))
        for pt in points:
            p.drawEllipse(pt, 3.0, 3.0)
