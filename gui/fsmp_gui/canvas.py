"""Interactive 2D canvas of the molecule editor.

Scene coordinates are angstroms in the xy plane; the view is mirrored
vertically so that +y points up. Atoms are draggable circles kept in sync
with the AtomTableModel (which is the single source of truth).
"""

import math
from enum import Enum

from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QBrush, QColor, QFont, QPainter, QPen, QTransform
from PySide6.QtWidgets import (QGraphicsEllipseItem, QGraphicsItem,
                               QGraphicsLineItem, QGraphicsScene,
                               QGraphicsSimpleTextItem, QGraphicsView)

from . import theme
from .atom_table import AtomTableModel
from .elements import BOND_TOLERANCE, covalent_radius, element_color
from .molecule import Atom


class Mode(Enum):
    SELECT = 0
    ADD = 1
    DELETE = 2


def display_radius(element: str) -> float:
    return min(max(0.45 * covalent_radius(element) + 0.18, 0.30), 0.85)


def text_color_for(background: str) -> QColor:
    c = QColor(background)
    lum = 0.299 * c.red() + 0.587 * c.green() + 0.114 * c.blue()
    return QColor("#10181a") if lum > 140 else QColor("#f2f5f7")


class AtomItem(QGraphicsEllipseItem):
    def __init__(self, canvas: "MoleculeCanvas", row: int):
        super().__init__()
        self.canvas = canvas
        self.row = row
        self.setFlags(QGraphicsItem.ItemIsMovable
                      | QGraphicsItem.ItemIsSelectable
                      | QGraphicsItem.ItemSendsScenePositionChanges)
        self.setZValue(1)
        self.label = QGraphicsSimpleTextItem(self)
        self.refresh()

    def refresh(self) -> None:
        atom = self.canvas.model.atom(self.row)
        r = display_radius(atom.element)
        color = QColor(element_color(atom.element))
        self.setRect(QRectF(-r, -r, 2 * r, 2 * r))
        self.setBrush(QBrush(color))
        self.setPen(QPen(color.darker(160), 0.05))
        self.setPos(atom.x, atom.y)

        self.label.setText(atom.element)
        font = QFont("Segoe UI", 10)
        font.setBold(True)
        self.label.setFont(font)
        self.label.setBrush(QBrush(text_color_for(element_color(atom.element))))
        bounds = self.label.boundingRect()
        k = 1.1 * r / max(bounds.height(), 1e-6)
        # mirror the text back (the view flips y) and center it on the atom
        self.label.setTransform(QTransform(k, 0, 0, -k, -k * bounds.width() / 2,
                                           k * bounds.height() / 2))

    def paint(self, painter, option, widget=None):
        painter.setRenderHint(QPainter.Antialiasing)
        painter.setBrush(self.brush())
        if self.isSelected():
            painter.setPen(QPen(QColor(theme.ACCENT_HOVER), 0.10))
        else:
            painter.setPen(self.pen())
        painter.drawEllipse(self.rect())

    def itemChange(self, change, value):
        if change == QGraphicsItem.ItemScenePositionHasChanged:
            self.canvas.atom_dragged(self)
        return super().itemChange(change, value)


class MoleculeCanvas(QGraphicsView):
    cursorMoved = Signal(float, float)
    atomAddRequested = Signal(float, float)  # click on empty space in ADD mode

    def __init__(self, model: AtomTableModel, parent=None):
        super().__init__(parent)
        self.model = model
        self.mode = Mode.SELECT
        self._syncing = False
        self._atom_items: list[AtomItem] = []
        self._bond_items: list[QGraphicsLineItem] = []

        self.setScene(QGraphicsScene(self))
        self.scene().setSceneRect(-60, -60, 120, 120)
        self.setRenderHint(QPainter.Antialiasing)
        self.setDragMode(QGraphicsView.RubberBandDrag)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.setMouseTracking(True)
        self.setBackgroundBrush(QColor(theme.BG_PANEL))
        self.setFrameShape(QGraphicsView.NoFrame)
        self.setTransform(QTransform.fromScale(46, -46))

        model.modelReset.connect(self.rebuild)
        model.rowsInserted.connect(self.rebuild)
        model.rowsRemoved.connect(self.rebuild)
        model.dataChanged.connect(self._model_edited)
        self.rebuild()

    # -- background grid ---------------------------------------------------

    def drawBackground(self, painter, rect):
        super().drawBackground(painter, rect)
        minor = QPen(QColor(theme.BORDER), 0)
        minor.setCosmetic(True)
        major = QPen(QColor(theme.BORDER).lighter(135), 0)
        major.setCosmetic(True)
        step = 1.0
        x0, x1 = math.floor(rect.left()), math.ceil(rect.right())
        y0, y1 = math.floor(rect.top()), math.ceil(rect.bottom())
        x = x0
        while x <= x1:
            painter.setPen(major if x % 5 == 0 else minor)
            painter.drawLine(QPointF(x, rect.top()), QPointF(x, rect.bottom()))
            x += step
        y = y0
        while y <= y1:
            painter.setPen(major if y % 5 == 0 else minor)
            painter.drawLine(QPointF(rect.left(), y), QPointF(rect.right(), y))
            y += step
        axes = QPen(QColor(theme.TEXT_DIM), 0)
        axes.setCosmetic(True)
        painter.setPen(axes)
        painter.drawLine(QPointF(rect.left(), 0), QPointF(rect.right(), 0))
        painter.drawLine(QPointF(0, rect.top()), QPointF(0, rect.bottom()))

    # -- model -> scene ----------------------------------------------------

    def rebuild(self) -> None:
        self._syncing = True
        self.scene().clear()  # also deletes label children
        self._atom_items = []
        self._bond_items = []
        for row in range(self.model.rowCount()):
            item = AtomItem(self, row)
            self.scene().addItem(item)
            self._atom_items.append(item)
        self._syncing = False
        self.rebuild_bonds()

    def _model_edited(self, top_left, bottom_right, roles=None) -> None:
        if self._syncing:
            return
        self._syncing = True
        for row in range(top_left.row(), bottom_right.row() + 1):
            self._atom_items[row].refresh()
        self._syncing = False
        self.rebuild_bonds()

    def rebuild_bonds(self) -> None:
        for line in self._bond_items:
            self.scene().removeItem(line)
        self._bond_items = []
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
                    self._bond_items.append(line)

    # -- scene -> model ----------------------------------------------------

    def atom_dragged(self, item: AtomItem) -> None:
        if self._syncing:
            return
        self._syncing = True
        pos = item.pos()
        self.model.move_atom(item.row, round(pos.x(), 4), round(pos.y(), 4))
        self._syncing = False
        self.rebuild_bonds()

    # -- interaction -------------------------------------------------------

    def set_mode(self, mode: Mode) -> None:
        self.mode = mode
        movable = mode == Mode.SELECT
        for item in self._atom_items:
            item.setFlag(QGraphicsItem.ItemIsMovable, movable)
        self.setDragMode(QGraphicsView.RubberBandDrag if movable
                         else QGraphicsView.NoDrag)

    def selected_rows(self) -> list[int]:
        return [it.row for it in self._atom_items if it.isSelected()]

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            pos = self.mapToScene(event.position().toPoint())
            hit = self.itemAt(event.position().toPoint())
            if isinstance(hit, QGraphicsSimpleTextItem):
                hit = hit.parentItem()
            if isinstance(hit, QGraphicsLineItem):
                hit = None  # bonds are not clickable
            if self.mode == Mode.ADD and hit is None:
                self.atomAddRequested.emit(round(pos.x(), 2), round(pos.y(), 2))
                return
            if self.mode == Mode.DELETE and isinstance(hit, AtomItem):
                self.model.remove_rows([hit.row])
                return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        pos = self.mapToScene(event.position().toPoint())
        self.cursorMoved.emit(pos.x(), pos.y())
        super().mouseMoveEvent(event)

    def wheelEvent(self, event):
        factor = 1.15 if event.angleDelta().y() > 0 else 1 / 1.15
        current = abs(self.transform().m11())
        if 4 < current * factor < 400:
            self.scale(factor, factor)

    def keyPressEvent(self, event):
        if event.key() in (Qt.Key_Delete, Qt.Key_Backspace):
            rows = self.selected_rows()
            if rows:
                self.model.remove_rows(rows)
                return
        super().keyPressEvent(event)

    def reset_view(self) -> None:
        atoms = self.model.molecule.atoms
        if atoms:
            xs = [a.x for a in atoms]
            ys = [a.y for a in atoms]
            rect = QRectF(QPointF(min(xs) - 2, min(ys) - 2),
                          QPointF(max(xs) + 2, max(ys) + 2))
        else:
            rect = QRectF(-6, -6, 12, 12)
        scale = min(self.viewport().width() / max(rect.width(), 1e-6),
                    self.viewport().height() / max(rect.height(), 1e-6))
        scale = min(max(scale, 4), 120)
        self.setTransform(QTransform.fromScale(scale, -scale))
        self.centerOn(rect.center())

    def add_atom_at(self, element: str, x: float, y: float) -> None:
        self.model.add_atom(Atom(element, x, y))
