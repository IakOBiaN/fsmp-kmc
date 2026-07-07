"""Qt table model over a SiteModel: shared source of truth for the site
editor's canvas and table."""

from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt, Signal

from .sitemodel import Site, SiteModel

COLUMNS = ["Label", "x, Å", "y, Å", "z, Å", "q, e", "ε, K", "σ, Å"]
_ATTR = ("label", "x", "y", "z", "q", "epsilon", "sigma")


class SiteTableModel(QAbstractTableModel):
    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._model = SiteModel()

    @property
    def site_model(self) -> SiteModel:
        return self._model

    def set_site_model(self, model: SiteModel) -> None:
        self.beginResetModel()
        self._model = model
        self.endResetModel()
        self.changed.emit()

    def site(self, row: int) -> Site:
        return self._model.sites[row]

    def add_site(self, site: Site) -> None:
        n = len(self._model.sites)
        self.beginInsertRows(QModelIndex(), n, n)
        self._model.sites.append(site)
        self.endInsertRows()
        self.changed.emit()

    def remove_rows(self, rows: list[int]) -> None:
        for row in sorted(set(rows), reverse=True):
            self.beginRemoveRows(QModelIndex(), row, row)
            del self._model.sites[row]
            self.endRemoveRows()
        if rows:
            self.changed.emit()

    def move_site(self, row: int, x: float, y: float) -> None:
        s = self._model.sites[row]
        s.x, s.y = x, y
        self.dataChanged.emit(self.index(row, 1), self.index(row, 2),
                              [Qt.DisplayRole, Qt.EditRole])
        self.changed.emit()

    # -- QAbstractTableModel ------------------------------------------------

    def rowCount(self, parent=QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._model.sites)

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
        s = self._model.sites[index.row()]
        col = index.column()
        if col == 0:
            return s.label
        value = getattr(s, _ATTR[col])
        if role == Qt.EditRole:
            return value
        if col in (4,):
            return f"{value:+.4f}"
        return f"{value:.4f}"

    def setData(self, index, value, role=Qt.EditRole):
        if role != Qt.EditRole:
            return False
        s = self._model.sites[index.row()]
        col = index.column()
        if col == 0:
            text = str(value).strip()
            if not text:
                return False
            s.label = text
        else:
            try:
                num = float(value)
            except (TypeError, ValueError):
                return False
            if col in (5, 6) and num < 0:  # epsilon, sigma are non-negative
                return False
            setattr(s, _ATTR[col], num)
        self.dataChanged.emit(index, index, [Qt.DisplayRole, Qt.EditRole])
        self.changed.emit()
        return True
