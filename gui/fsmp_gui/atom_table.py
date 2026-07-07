"""Qt table model over a Molecule: the single source of truth that the
canvas and the table view both edit and listen to."""

from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt, Signal

from .elements import normalize_symbol
from .molecule import Atom, Molecule

COLUMNS = ["Element", "x, Å", "y, Å", "z, Å"]


class AtomTableModel(QAbstractTableModel):
    changed = Signal()  # any modification, structural or not

    def __init__(self, parent=None):
        super().__init__(parent)
        self._mol = Molecule()

    # -- molecule access ---------------------------------------------------

    @property
    def molecule(self) -> Molecule:
        return self._mol

    def set_molecule(self, mol: Molecule) -> None:
        self.beginResetModel()
        self._mol = mol
        self.endResetModel()
        self.changed.emit()

    def atom(self, row: int) -> Atom:
        return self._mol.atoms[row]

    def add_atom(self, atom: Atom) -> None:
        n = len(self._mol.atoms)
        self.beginInsertRows(QModelIndex(), n, n)
        self._mol.atoms.append(atom)
        self.endInsertRows()
        self.changed.emit()

    def remove_rows(self, rows: list[int]) -> None:
        for row in sorted(set(rows), reverse=True):
            self.beginRemoveRows(QModelIndex(), row, row)
            del self._mol.atoms[row]
            self.endRemoveRows()
        if rows:
            self.changed.emit()

    def move_atom(self, row: int, x: float, y: float) -> None:
        """Position update coming from a canvas drag."""
        a = self._mol.atoms[row]
        a.x, a.y = x, y
        tl, br = self.index(row, 1), self.index(row, 2)
        self.dataChanged.emit(tl, br, [Qt.DisplayRole, Qt.EditRole])
        self.changed.emit()

    # -- QAbstractTableModel -----------------------------------------------

    def rowCount(self, parent=QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._mol.atoms)

    def columnCount(self, parent=QModelIndex()) -> int:
        return len(COLUMNS)

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return None
        if orientation == Qt.Horizontal:
            return COLUMNS[section]
        return section + 1

    def flags(self, index):
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable

    def data(self, index, role=Qt.DisplayRole):
        if role not in (Qt.DisplayRole, Qt.EditRole):
            return None
        a = self._mol.atoms[index.row()]
        col = index.column()
        if col == 0:
            return a.element
        value = (a.x, a.y, a.z)[col - 1]
        return value if role == Qt.EditRole else f"{value:.4f}"

    def setData(self, index, value, role=Qt.EditRole):
        if role != Qt.EditRole:
            return False
        a = self._mol.atoms[index.row()]
        col = index.column()
        if col == 0:
            symbol = normalize_symbol(str(value))
            if not symbol or not symbol[0].isalpha():
                return False
            a.element = symbol
        else:
            try:
                num = float(value)
            except (TypeError, ValueError):
                return False
            if col == 1:
                a.x = num
            elif col == 2:
                a.y = num
            else:
                a.z = num
        self.dataChanged.emit(index, index, [Qt.DisplayRole, Qt.EditRole])
        self.changed.emit()
        return True
