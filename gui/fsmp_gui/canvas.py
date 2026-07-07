"""Interactive 2D editor canvases.

ItemCanvas is the shared machinery (draggable circular items synced to a
table model, add/select/delete modes); MoleculeCanvas and SiteCanvas fill in
what an item looks like and where it lives. Scene coordinates are angstroms
with +y up (the view flips y, so item labels are mirrored back).
"""

import math
from enum import Enum

from PySide6.QtCore import QRectF, Qt, Signal
from PySide6.QtGui import QBrush, QColor, QFont, QPainter, QPen, QTransform
from PySide6.QtWidgets import (QGraphicsEllipseItem, QGraphicsItem,
                               QGraphicsLineItem, QGraphicsSimpleTextItem)

from . import theme
from .atom_table import AtomTableModel
from .elements import BOND_TOLERANCE, covalent_radius, element_color
from .gridview import GridView
from .molecule import Atom
from .site_table import SiteTableModel


class Mode(Enum):
    SELECT = 0
    ADD = 1
    DELETE = 2


def _text_color_for(color: QColor) -> QColor:
    lum = 0.299 * color.red() + 0.587 * color.green() + 0.114 * color.blue()
    return QColor("#10181a") if lum > 140 else QColor("#f2f5f7")


class CircleItem(QGraphicsEllipseItem):
    """A draggable filled circle with a centered, view-corrected label.
    Subclasses supply radius/fill/text/ring from the model row."""

    def __init__(self, canvas: "ItemCanvas", row: int):
        super().__init__()
        self.canvas = canvas
        self.row = row
        flags = (QGraphicsItem.ItemIsSelectable
                 | QGraphicsItem.ItemSendsScenePositionChanges)
        if not canvas.read_only:
            flags |= QGraphicsItem.ItemIsMovable
        self.setFlags(flags)
        self.setZValue(1)
        self.label = QGraphicsSimpleTextItem(self)
        self.refresh()

    # subclass hooks
    def radius(self) -> float: ...
    def fill(self) -> QColor: ...
    def text(self) -> str: ...
    def ring(self) -> QColor | None:
        return None

    def refresh(self) -> None:
        r = self.radius()
        color = self.fill()
        self.setRect(QRectF(-r, -r, 2 * r, 2 * r))
        self.setBrush(QBrush(color))
        self.setPos(*self.canvas.pos(self.row))
        self.label.setText(self.text())
        font = QFont("Segoe UI", 10)
        font.setBold(True)
        self.label.setFont(font)
        self.label.setBrush(QBrush(_text_color_for(color)))
        b = self.label.boundingRect()
        k = 1.0 * r / max(b.height(), 1e-6)
        self.label.setTransform(QTransform(k, 0, 0, -k, -k * b.width() / 2,
                                           k * b.height() / 2))

    def paint(self, painter, option, widget=None):
        painter.setRenderHint(QPainter.Antialiasing)
        painter.setBrush(self.brush())
        if self.isSelected():
            painter.setPen(QPen(QColor(theme.ACCENT_HOVER), 0.10))
        elif (ring := self.ring()) is not None:
            painter.setPen(QPen(ring, 0.07))
        else:
            painter.setPen(QPen(self.brush().color().darker(160), 0.05))
        painter.drawEllipse(self.rect())

    def itemChange(self, change, value):
        if change == QGraphicsItem.ItemScenePositionHasChanged:
            self.canvas.item_dragged(self)
        return super().itemChange(change, value)


class ItemCanvas(GridView):
    addRequested = Signal(float, float)  # empty-space click in ADD mode

    def __init__(self, model, parent=None, read_only=False):
        super().__init__(parent)
        self.model = model
        self.mode = Mode.SELECT
        self.read_only = read_only
        self._syncing = False
        self._items: list[CircleItem] = []
        model.modelReset.connect(self.rebuild)
        model.rowsInserted.connect(self.rebuild)
        model.rowsRemoved.connect(self.rebuild)
        model.dataChanged.connect(self._model_edited)
        self.rebuild()

    # subclass hooks
    def make_item(self, row: int) -> CircleItem: ...
    def pos(self, row: int) -> tuple[float, float]: ...
    def move(self, row: int, x: float, y: float) -> None: ...
    def points(self) -> list[tuple[float, float]]: ...
    def decorate(self) -> None:
        """Extra scene items (e.g. bonds); rebuilt on any change."""

    def rebuild(self) -> None:
        self._syncing = True
        self.scene().clear()
        self._items = [self.make_item(row) for row in range(self.model.rowCount())]
        for item in self._items:
            self.scene().addItem(item)
        self._syncing = False
        self.decorate()

    def _model_edited(self, top_left, bottom_right, roles=None) -> None:
        if self._syncing:
            return
        self._syncing = True
        for row in range(top_left.row(), bottom_right.row() + 1):
            self._items[row].refresh()
        self._syncing = False
        self.decorate()

    def item_dragged(self, item: CircleItem) -> None:
        if self._syncing:
            return
        self._syncing = True
        p = item.pos()
        self.move(item.row, round(p.x(), 4), round(p.y(), 4))
        self._syncing = False
        self.decorate()

    def set_mode(self, mode: Mode) -> None:
        if self.read_only:
            return
        self.mode = mode
        movable = mode == Mode.SELECT
        for item in self._items:
            item.setFlag(QGraphicsItem.ItemIsMovable, movable)
        self.setDragMode(GridView.RubberBandDrag if movable else GridView.NoDrag)

    def selected_rows(self) -> list[int]:
        return [it.row for it in self._items if it.isSelected()]

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            pos = self.mapToScene(event.position().toPoint())
            hit = self.itemAt(event.position().toPoint())
            if isinstance(hit, QGraphicsSimpleTextItem):
                hit = hit.parentItem()
            if not isinstance(hit, CircleItem):
                hit = None  # empty space, grid or a bond
            if self.mode == Mode.ADD and hit is None:
                self.addRequested.emit(round(pos.x(), 2), round(pos.y(), 2))
                return
            if self.mode == Mode.DELETE and hit is not None:
                self.model.remove_rows([hit.row])
                return
        super().mousePressEvent(event)

    def keyPressEvent(self, event):
        if not self.read_only and event.key() in (Qt.Key_Delete, Qt.Key_Backspace):
            if self.selected_rows():
                self.model.remove_rows(self.selected_rows())
                return
        super().keyPressEvent(event)

    def reset_view(self) -> None:
        self.fit_points(self.points())


# -- atomistic molecule ----------------------------------------------------

def _atom_radius(element: str) -> float:
    return min(max(0.45 * covalent_radius(element) + 0.18, 0.30), 0.85)


class AtomItem(CircleItem):
    def radius(self):
        return _atom_radius(self.canvas.model.atom(self.row).element)

    def fill(self):
        return QColor(element_color(self.canvas.model.atom(self.row).element))

    def text(self):
        return self.canvas.model.atom(self.row).element


class MoleculeCanvas(ItemCanvas):
    def __init__(self, model: AtomTableModel, parent=None, read_only=False):
        self._bonds: list[QGraphicsLineItem] = []
        super().__init__(model, parent, read_only)

    def make_item(self, row):
        return AtomItem(self, row)

    def pos(self, row):
        a = self.model.atom(row)
        return a.x, a.y

    def move(self, row, x, y):
        self.model.move_atom(row, x, y)

    def points(self):
        return [(a.x, a.y) for a in self.model.molecule.atoms]

    def decorate(self):
        for line in self._bonds:
            self.scene().removeItem(line)
        self._bonds = []
        atoms = self.model.molecule.atoms
        pen = QPen(QColor("#5b6f70"), 0.16, Qt.SolidLine, Qt.RoundCap)
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                a, b = atoms[i], atoms[j]
                d = math.dist((a.x, a.y, a.z), (b.x, b.y, b.z))
                limit = BOND_TOLERANCE * (covalent_radius(a.element)
                                          + covalent_radius(b.element))
                if 1e-6 < d < limit:
                    line = QGraphicsLineItem(a.x, a.y, b.x, b.y)
                    line.setPen(pen)
                    line.setZValue(0)
                    self.scene().addItem(line)
                    self._bonds.append(line)

    def add_atom_at(self, element: str, x: float, y: float) -> None:
        self.model.add_atom(Atom(element, x, y))


# -- site model ------------------------------------------------------------

def _site_radius(sigma: float) -> float:
    return 0.34 if sigma <= 0 else min(max(0.22 * sigma + 0.12, 0.32), 0.95)


def site_color(q: float) -> QColor:
    return QColor("#4a78f0" if q > 0 else "#f0554a" if q < 0 else "#6f8a86")


class SiteItem(CircleItem):
    def radius(self):
        return _site_radius(self.canvas.model.site(self.row).sigma)

    def fill(self):
        return site_color(self.canvas.model.site(self.row).q)

    def text(self):
        return self.canvas.model.site(self.row).label

    def ring(self):
        return QColor(theme.GOLD) if self.canvas.model.site(self.row).is_lj else None


class SiteCanvas(ItemCanvas):
    def __init__(self, model: SiteTableModel, parent=None, read_only=False):
        super().__init__(model, parent, read_only)

    def make_item(self, row):
        return SiteItem(self, row)

    def pos(self, row):
        s = self.model.site(row)
        return s.x, s.y

    def move(self, row, x, y):
        self.model.move_site(row, x, y)

    def points(self):
        return [(s.x, s.y) for s in self.model.site_model.sites]
