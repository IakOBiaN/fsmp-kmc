"""A molecule placement inside the unit cell (position and orientation) and a
Qt table model over the list of placements, shared by the cell canvas and the
placements table."""

from dataclasses import dataclass

from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt, Signal

COLUMNS = ["x, Å", "y, Å", "φ, °"]
_ATTR = ("x", "y", "phi")


@dataclass
class Placement:
    x: float
    y: float
    phi: float = 0.0

    def to_dict(self) -> dict:
        return {"x": self.x, "y": self.y, "phi": self.phi}

    @staticmethod
    def from_dict(d: dict) -> "Placement":
        return Placement(float(d["x"]), float(d["y"]), float(d.get("phi", 0.0)))


class PlacementTableModel(QAbstractTableModel):
    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._items: list[Placement] = []

    @property
    def placements(self) -> list[Placement]:
        return self._items

    def set_placements(self, items: list[Placement]) -> None:
        self.beginResetModel()
        self._items = list(items)
        self.endResetModel()
        self.changed.emit()

    def at(self, row: int) -> Placement:
        return self._items[row]

    def add(self, p: Placement) -> None:
        n = len(self._items)
        self.beginInsertRows(QModelIndex(), n, n)
        self._items.append(p)
        self.endInsertRows()
        self.changed.emit()

    def remove_rows(self, rows: list[int]) -> None:
        for row in sorted(set(rows), reverse=True):
            if 0 <= row < len(self._items):
                self.beginRemoveRows(QModelIndex(), row, row)
                del self._items[row]
                self.endRemoveRows()
        if rows:
            self.changed.emit()

    def move(self, row: int, x: float, y: float) -> None:
        p = self._items[row]
        p.x, p.y = x, y
        self.dataChanged.emit(self.index(row, 0), self.index(row, 1),
                              [Qt.DisplayRole, Qt.EditRole])
        self.changed.emit()

    def set_phi(self, row: int, phi_deg: float) -> None:
        """Orientation update coming from the canvas rotation handle; keeps
        phi in [0, 360)."""
        p = self._items[row]
        p.phi = phi_deg % 360.0
        idx = self.index(row, 2)
        self.dataChanged.emit(idx, idx, [Qt.DisplayRole, Qt.EditRole])
        self.changed.emit()

    # -- QAbstractTableModel ------------------------------------------------

    def rowCount(self, parent=QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._items)

    def columnCount(self, parent=QModelIndex()) -> int:
        return len(COLUMNS)

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return None
        return COLUMNS[section] if orientation == Qt.Horizontal else section + 1

    def flags(self, index):
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable

    def data(self, index, role=Qt.DisplayRole):
        if role not in (Qt.DisplayRole, Qt.EditRole):
            return None
        value = getattr(self._items[index.row()], _ATTR[index.column()])
        return value if role == Qt.EditRole else f"{value:.4f}"

    def setData(self, index, value, role=Qt.EditRole):
        if role != Qt.EditRole:
            return False
        try:
            num = float(value)
        except (TypeError, ValueError):
            return False
        setattr(self._items[index.row()], _ATTR[index.column()], num)
        self.dataChanged.emit(index, index, [Qt.DisplayRole, Qt.EditRole])
        self.changed.emit()
        return True
