"""Tab 4 "Unit cell": lay out molecule copies in a rectangular cell, resize
the cell and move the copies. Optimization (via the engine) comes later; this
tab only builds and stores the rough cell.

Every copy is an instance of the project molecule model; the cell stores each
copy's position and orientation.
"""

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QAbstractItemView, QButtonGroup, QDoubleSpinBox,
                               QFrame, QHBoxLayout, QHeaderView, QLabel,
                               QMessageBox, QPushButton, QSplitter, QTableView,
                               QToolButton, QVBoxLayout, QWidget)

from ..canvas import Mode
from ..glyph import model_glyph
from ..placement_table import Placement, PlacementTableModel
from ..project import Project
from ..unit_cell_canvas import UnitCellCanvas


class UnitCellTab(QWidget):
    statusMessage = Signal(str)

    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self.project = project
        self.model = PlacementTableModel(self)
        self.model.changed.connect(self._on_changed)
        self._dirty = False

        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.addLayout(self._build_toolbar())

        splitter = QSplitter(Qt.Horizontal)
        root.addWidget(splitter, 1)
        self.canvas = UnitCellCanvas(self.model)
        self.canvas.cursorMoved.connect(
            lambda x, y: self.statusMessage.emit(f"x = {x:.2f} Å,  y = {y:.2f} Å"))
        splitter.addWidget(self.canvas)
        splitter.addWidget(self._build_side())
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 0)
        splitter.setSizes([880, 380])

        self._set_mode(Mode.SELECT)

    # -- construction ------------------------------------------------------

    def _build_toolbar(self) -> QHBoxLayout:
        bar = QHBoxLayout()
        bar.setSpacing(6)

        self.mode_group = QButtonGroup(self)
        self.mode_buttons = {}
        for mode, text, tip in ((Mode.SELECT, "Select", "Drag molecules"),
                                (Mode.ADD, "Add", "Click to add a molecule copy"),
                                (Mode.DELETE, "Delete", "Click a molecule to remove it")):
            b = QToolButton()
            b.setText(text)
            b.setToolTip(tip)
            b.setCheckable(True)
            b.clicked.connect(lambda _=False, m=mode: self._set_mode(m))
            self.mode_group.addButton(b)
            self.mode_buttons[mode] = b
            bar.addWidget(b)

        bar.addSpacing(16)
        bar.addWidget(QLabel("cell x"))
        self.cell_x = self._cell_spin(12.0)
        bar.addWidget(self.cell_x)
        bar.addWidget(QLabel("cell y"))
        self.cell_y = self._cell_spin(12.0)
        bar.addWidget(self.cell_y)

        bar.addSpacing(12)
        fit = QPushButton("Fit view")
        fit.clicked.connect(lambda: self.canvas._redraw_fit())
        bar.addWidget(fit)
        bar.addStretch(1)

        attach = QPushButton("Use in project")
        attach.setProperty("primary", True)
        attach.clicked.connect(self.use_in_project)
        bar.addWidget(attach)
        return bar

    def _cell_spin(self, value: float) -> QDoubleSpinBox:
        s = QDoubleSpinBox()
        s.setRange(1.0, 200.0)
        s.setDecimals(3)
        s.setSingleStep(0.5)
        s.setSuffix(" Å")
        s.setValue(value)
        s.valueChanged.connect(self._cell_changed)
        return s

    def _build_side(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(8, 0, 0, 0)

        self.info = QLabel()
        self.info.setWordWrap(True)
        layout.addWidget(self.info)

        caption = QLabel("Molecules in the cell")
        caption.setProperty("dim", True)
        layout.addWidget(caption)

        self.table = QTableView()
        self.table.setModel(self.model)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setAlternatingRowColors(True)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.verticalHeader().setDefaultSectionSize(26)
        layout.addWidget(self.table, 1)

        row = QHBoxLayout()
        add = QPushButton("Add molecule")
        add.clicked.connect(self._add_center)
        rem = QPushButton("Remove selected")
        rem.clicked.connect(self._remove_selected)
        row.addWidget(add)
        row.addWidget(rem)
        row.addStretch(1)
        layout.addLayout(row)

        layout.addSpacing(8)
        caption2 = QLabel("Project unit cell")
        caption2.setProperty("dim", True)
        layout.addWidget(caption2)
        frame = QFrame()
        frame.setProperty("panel", True)
        fl = QVBoxLayout(frame)
        self.slot_label = QLabel()
        self.slot_label.setWordWrap(True)
        fl.addWidget(self.slot_label)
        srow = QHBoxLayout()
        self.slot_load = QPushButton("Load into editor")
        self.slot_load.clicked.connect(self._load_from_project)
        self.slot_clear = QPushButton("Remove")
        self.slot_clear.clicked.connect(self._clear_project)
        srow.addWidget(self.slot_load)
        srow.addWidget(self.slot_clear)
        srow.addStretch(1)
        fl.addLayout(srow)
        layout.addWidget(frame)
        return panel

    # -- state -------------------------------------------------------------

    def refresh(self) -> None:
        """Called when the tab becomes current: pick up the project model glyph
        and load the stored cell (or start a default one)."""
        self.canvas.set_glyph(model_glyph(self.project))
        if not self.model.placements and not self._dirty:
            if self.project.unit_cell is not None:
                self._load_from_project(confirm=False)
            else:
                self._start_default()
        self._refresh_slot()
        self._update_info()

    def _start_default(self) -> None:
        self.cell_x.setValue(12.0)
        self.cell_y.setValue(12.0)
        self.canvas.set_cell(12.0, 12.0)
        self.model.set_placements([Placement(6.0, 6.0, 0.0)])
        self._dirty = False

    def _on_changed(self) -> None:
        self._dirty = True
        self._update_info()

    def _cell_changed(self) -> None:
        self.canvas.set_cell(self.cell_x.value(), self.cell_y.value())
        self._dirty = True
        self._update_info()

    def _update_info(self) -> None:
        n = len(self.model.placements)
        area = self.cell_x.value() * self.cell_y.value()
        density = (n / area * 100.0) if area else 0.0  # molecules per nm^2
        self.info.setText(f"{n} molecules   ·   cell {self.cell_x.value():.2f} × "
                          f"{self.cell_y.value():.2f} Å   ·   {density:.3f} nm⁻²")

    def _set_mode(self, mode: Mode) -> None:
        self.mode_buttons[mode].setChecked(True)
        self.canvas.set_mode(mode)

    # -- editing -----------------------------------------------------------

    def _add_center(self) -> None:
        self.model.add(Placement(round(self.cell_x.value() / 2, 3),
                                 round(self.cell_y.value() / 2, 3), 0.0))

    def _remove_selected(self) -> None:
        rows = [i.row() for i in self.table.selectionModel().selectedRows()]
        self.model.remove_rows(rows)

    # -- project integration -----------------------------------------------

    def _refresh_slot(self) -> None:
        uc = self.project.unit_cell
        if uc is None:
            self.slot_label.setText("Not saved yet. Build the cell and press "
                                    "'Use in project'.")
            self.slot_load.setEnabled(False)
            self.slot_clear.setEnabled(False)
        else:
            self.slot_label.setText(f"{len(uc['molecules'])} molecules, cell "
                                    f"{uc['cell_x']:.2f} × {uc['cell_y']:.2f} Å")
            self.slot_load.setEnabled(True)
            self.slot_clear.setEnabled(True)

    def use_in_project(self) -> None:
        if not self.model.placements:
            QMessageBox.information(self, "Empty cell", "Add at least one molecule.")
            return
        self.project.set_unit_cell(self.cell_x.value(), self.cell_y.value(),
                                   [p.to_dict() for p in self.model.placements])
        self._dirty = False
        self._refresh_slot()
        self.statusMessage.emit("Unit cell saved to the project")

    def _load_from_project(self, confirm: bool = True) -> None:
        uc = self.project.unit_cell
        if uc is None:
            return
        if confirm and self._dirty and self.model.placements:
            if QMessageBox.question(self, "Discard changes",
                                    "Replace the current cell with the saved one?"
                                    ) != QMessageBox.Yes:
                return
        self.cell_x.setValue(uc["cell_x"])
        self.cell_y.setValue(uc["cell_y"])
        self.canvas.set_cell(uc["cell_x"], uc["cell_y"])
        self.model.set_placements([Placement.from_dict(d) for d in uc["molecules"]])
        self._dirty = False
        self._update_info()

    def _clear_project(self) -> None:
        if self.project.unit_cell is None:
            return
        if QMessageBox.question(self, "Remove unit cell",
                                "Remove the saved unit cell from the project?"
                                ) != QMessageBox.Yes:
            return
        self.project.clear_unit_cell()
        self._refresh_slot()
        self.statusMessage.emit("Unit cell removed from the project")
