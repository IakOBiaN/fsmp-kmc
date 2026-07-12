"""Lazy viewer for the run trajectory (extended xyz, possibly hundreds of
megabytes). It always shows the newest complete frame: a background thread
indexes frame byte offsets incrementally, poll() picks up whatever the
engine has appended since, and only the current frame is ever read into
memory, so the interface stays light no matter the file size.

Atoms are drawn as element-colored dots; the opacity follows the lambda
column, so gas-phase (damped) molecules fade out exactly where the field
says they should."""

import os
import re

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtGui import QBrush, QColor, QPen
from PySide6.QtWidgets import (QGraphicsRectItem, QGraphicsView, QHBoxLayout,
                               QLabel, QPushButton, QVBoxLayout, QWidget)

from . import theme
from .elements import covalent_radius, element_color
from .gridview import GridView

_LATTICE = re.compile(r'Lattice="([^"]+)"')


class FrameIndexer(QThread):
    """Collect the byte offset of every complete frame. Incremental: pass
    the previous offsets and end position to continue after the file has
    grown. A frame the engine is still writing is left for the next pass."""

    indexed = Signal(list, int)   # offsets, end position

    def __init__(self, path, offsets=None, pos=0, parent=None):
        super().__init__(parent)
        self._path = path
        self._offsets = list(offsets or [])
        self._pos = pos
        self._cancel = False

    def cancel(self):
        self._cancel = True

    def run(self):
        offsets, pos = self._offsets, self._pos
        try:
            with open(self._path, "rb") as f:
                f.seek(pos)
                while not self._cancel:
                    start = f.tell()
                    head = f.readline()
                    if not head.endswith(b"\n"):
                        break
                    try:
                        count = int(head.split()[0])
                    except (ValueError, IndexError):
                        break
                    torn = not f.readline().endswith(b"\n")
                    for _ in range(count):
                        if torn:
                            break
                        torn = not f.readline().endswith(b"\n")
                    if torn:
                        break
                    offsets.append(start)
                    pos = f.tell()
        except OSError:
            pass
        self.indexed.emit(offsets, pos)


def read_frame(path, offset: int):
    """((lx, ly), [(element, x, y, lambda), ...]) of one frame."""
    with open(path, "rb") as f:
        f.seek(offset)
        count = int(f.readline().split()[0])
        comment = f.readline().decode("ascii", "replace")
        atoms = []
        for _ in range(count):
            cols = f.readline().split()
            lam = float(cols[4]) if len(cols) > 4 else 1.0
            atoms.append((cols[0].decode("ascii", "replace"),
                          float(cols[1]), float(cols[2]), lam))
    m = _LATTICE.search(comment)
    lattice = None
    if m:
        lat = [float(t) for t in m.group(1).split()]
        lattice = (lat[0], lat[4])
    return lattice, atoms


class TrajectoryCanvas(GridView):
    MIN_SCALE = 0.2
    SHOW_GRID = False

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setDragMode(QGraphicsView.ScrollHandDrag)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self._fitted = False

    def show_frame(self, lattice, atoms) -> None:
        scene = self.scene()
        scene.clear()
        if lattice:
            box = QGraphicsRectItem(0, 0, lattice[0], lattice[1])
            box.setPen(QPen(QColor(theme.ACCENT), 0))
            box.setZValue(-1)
            scene.addItem(box)
        pen = QPen(Qt.NoPen)
        for element, x, y, lam in atoms:
            r = min(max(0.4 * covalent_radius(element) + 0.15, 0.3), 0.8)
            dot = scene.addEllipse(x - r, y - r, 2 * r, 2 * r, pen,
                                   QBrush(QColor(element_color(element))))
            dot.setOpacity(0.25 + 0.75 * min(max(lam, 0.0), 1.0))
        if not self._fitted:
            self.fit_frame(lattice, atoms)

    def fit_frame(self, lattice, atoms) -> None:
        if lattice:
            points = [(0.0, 0.0), lattice]
        else:
            points = [(x, y) for _, x, y, _ in atoms] or [(0.0, 0.0)]
        self.fit_points(points, pad=4)
        self._fitted = True


class TrajectoryViewer(QWidget):
    """The newest frame of one trajectory file, kept fresh by poll()."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._path = None
        self._offsets: list = []
        self._end = 0
        self._shown = 0            # frame count at the last render
        self._indexer: FrameIndexer | None = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.canvas = TrajectoryCanvas()
        layout.addWidget(self.canvas, 1)

        bar = QHBoxLayout()
        self.info = QLabel("no trajectory yet")
        self.info.setProperty("dim", True)
        bar.addWidget(self.info)
        bar.addStretch(1)
        self.fit_btn = QPushButton("Fit")
        self.fit_btn.clicked.connect(self._refit)
        bar.addWidget(self.fit_btn)
        layout.addLayout(bar)

    def set_path(self, path) -> None:
        """Point the viewer at a trajectory file (resets the index)."""
        if path == self._path:
            return
        self._path = path
        self._offsets, self._end, self._shown = [], 0, 0
        self.canvas.scene().clear()
        self.canvas._fitted = False
        self.info.setText("no frames yet")
        self.poll()

    def poll(self) -> None:
        """Index frames appended since the last look and show the newest."""
        if self._path is None or self._indexer is not None:
            return
        try:
            size = os.path.getsize(self._path)
        except OSError:
            return
        if size <= self._end and self._offsets:
            return
        self._indexer = FrameIndexer(self._path, self._offsets, self._end,
                                     self)
        self._indexer.indexed.connect(self._on_indexed)
        self._indexer.finished.connect(lambda: setattr(self, "_indexer", None))
        self._indexer.start()

    def _on_indexed(self, offsets: list, end: int) -> None:
        self._offsets, self._end = offsets, end
        if not offsets:
            self.info.setText("no frames yet")
        elif len(offsets) != self._shown:
            self._show_last()

    def _show_last(self) -> None:
        try:
            lattice, atoms = read_frame(self._path, self._offsets[-1])
        except (OSError, ValueError, IndexError):
            return
        self.canvas.show_frame(lattice, atoms)
        self._shown = len(self._offsets)
        self.info.setText(f"frame {self._shown} (newest)")

    def _refit(self) -> None:
        self.canvas._fitted = False
        self._shown = 0
        if self._offsets:
            self._show_last()
