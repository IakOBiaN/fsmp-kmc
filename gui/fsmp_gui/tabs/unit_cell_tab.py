"""Tab 4 "Unit cell": lay out molecule copies in a rectangular cell, resize
the cell, move and rotate the copies (the canvas and the placements table
share one selection). The Optimize button hands the rough cell to
the engine (structure = calculate, optimize_only) and plays its animation
back live until the optimized cell lands in the editor. Cells are exchanged
with .cell files; the repository cells/ folder holds the reference
structures from StructureGenerator.h in that format.

Every copy is an instance of the project molecule model; the cell stores each
copy's position and orientation.
"""

from pathlib import Path

from PySide6.QtCore import QSettings, Qt, Signal
from PySide6.QtWidgets import (QAbstractItemView, QButtonGroup, QDoubleSpinBox,
                               QFileDialog, QFrame, QHBoxLayout, QHeaderView,
                               QLabel, QMessageBox, QPushButton, QSplitter,
                               QTableView, QToolButton, QVBoxLayout, QWidget)

from .. import cellfile
from ..canvas import Mode
from ..engine import EngineError, app_root, prepare_run
from ..glyph import model_glyph
from ..placement_table import Placement, PlacementTableModel
from ..project import Project
from ..unit_cell_canvas import UnitCellCanvas


class UnitCellTab(QWidget):
    statusMessage = Signal(str)
    projectCellChanged = Signal()

    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self.project = project
        self.model = PlacementTableModel(self)
        self.model.changed.connect(self._on_changed)
        self._dirty = False
        self._run = None

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

        # canvas and table always point at the same molecule
        self.canvas.moleculePicked.connect(self._on_canvas_pick)
        self.table.selectionModel().selectionChanged.connect(
            self._on_table_selection)

        self._set_mode(Mode.SELECT)

    # -- construction ------------------------------------------------------

    def _build_toolbar(self) -> QHBoxLayout:
        bar = QHBoxLayout()
        bar.setSpacing(6)

        for text, slot in (("Open…", self.open_cell),
                           ("Save as…", self.save_cell)):
            b = QPushButton(text)
            b.clicked.connect(slot)
            bar.addWidget(b)
        bar.addSpacing(16)

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

        bar.addSpacing(12)
        self.optimize_btn = QPushButton("Optimize")
        self.optimize_btn.setToolTip("Optimize the cell with the engine "
                                     "(uses the project potential)")
        self.optimize_btn.clicked.connect(self._optimize)
        bar.addWidget(self.optimize_btn)
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self._stop_optimize)
        bar.addWidget(self.stop_btn)
        bar.addStretch(1)

        self.attach_btn = QPushButton("Use in project")
        self.attach_btn.setProperty("primary", True)
        self.attach_btn.clicked.connect(self.use_in_project)
        bar.addWidget(self.attach_btn)
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

        self.opt_status = QLabel(" ")
        self.opt_status.setWordWrap(True)
        self.opt_status.setProperty("dim", True)
        layout.addWidget(self.opt_status)

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
        self.add_btn = QPushButton("Add molecule")
        self.add_btn.clicked.connect(self._add_center)
        self.rem_btn = QPushButton("Remove selected")
        self.rem_btn.clicked.connect(self._remove_selected)
        row.addWidget(self.add_btn)
        row.addWidget(self.rem_btn)
        row.addStretch(1)
        layout.addLayout(row)

        hint = QLabel("Click a molecule to select it, then drag its round "
                      "handle to rotate (hold Shift to snap to 5°); the φ "
                      "column takes exact values.")
        hint.setWordWrap(True)
        hint.setProperty("dim", True)
        layout.addWidget(hint)

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
        area = self.cell_x.value() * self.cell_y.value()   # A^2
        # the engine's density unit: micromoles per square meter
        density = (n / area * 1.0e26 / 6.02214076e23) if area else 0.0
        self.info.setText(f"{n} molecules   ·   cell {self.cell_x.value():.2f} × "
                          f"{self.cell_y.value():.2f} Å   ·   "
                          f"density {density:.4f} µmol/m²")

    def _set_mode(self, mode: Mode) -> None:
        self.mode_buttons[mode].setChecked(True)
        self.canvas.set_mode(mode)

    # -- selection sync ------------------------------------------------------

    def _on_canvas_pick(self, row: int) -> None:
        if row < 0:
            self.table.clearSelection()
        else:
            self.table.selectRow(row)

    def _on_table_selection(self, *_) -> None:
        rows = [i.row() for i in self.table.selectionModel().selectedRows()]
        self.canvas.set_selected(rows[0] if len(rows) == 1 else None)

    # -- cell files ----------------------------------------------------------

    def _cells_dir(self) -> str:
        stored = QSettings().value("cells/last_dir", "")
        if stored:
            return stored
        bundled = app_root() / "cells"
        return str(bundled if bundled.is_dir() else Path.home())

    def open_cell(self) -> None:
        if self._dirty and self.model.placements:
            if QMessageBox.question(self, "Discard changes",
                                    "Replace the current cell with the file?"
                                    ) != QMessageBox.Yes:
                return
        path, _ = QFileDialog.getOpenFileName(
            self, "Open unit cell", self._cells_dir(),
            f"Unit cells (*{cellfile.SUFFIX});;All files (*)")
        if not path:
            return
        try:
            cell_x, cell_y, placements, _ = cellfile.load_cell(path)
        except (OSError, ValueError) as e:
            QMessageBox.warning(self, "Cannot open unit cell", f"{path}\n\n{e}")
            return
        QSettings().setValue("cells/last_dir", str(Path(path).parent))
        self._apply_cell(cell_x, cell_y, placements)
        self._dirty = True   # differs from the project slot until saved there
        self._update_info()
        self.statusMessage.emit(f"Opened {path}")

    def save_cell(self) -> None:
        if not self.model.placements:
            QMessageBox.information(self, "Empty cell",
                                    "Add at least one molecule.")
            return
        suggestion = str(Path(self._cells_dir()) / f"unit_cell{cellfile.SUFFIX}")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save unit cell", suggestion,
            f"Unit cells (*{cellfile.SUFFIX})")
        if not path:
            return
        try:
            cellfile.save_cell(path, self.cell_x.value(), self.cell_y.value(),
                               [(p.x, p.y, p.phi) for p in self.model.placements])
        except OSError as e:
            QMessageBox.warning(self, "Cannot save unit cell", str(e))
            return
        QSettings().setValue("cells/last_dir", str(Path(path).parent))
        self.statusMessage.emit(f"Saved {path}")

    # -- editing -----------------------------------------------------------

    def _add_center(self) -> None:
        self.model.add(Placement(round(self.cell_x.value() / 2, 3),
                                 round(self.cell_y.value() / 2, 3), 0.0))

    def _remove_selected(self) -> None:
        rows = [i.row() for i in self.table.selectionModel().selectedRows()]
        self.model.remove_rows(rows)

    # -- optimization via the engine -----------------------------------------

    def _optimize(self) -> None:
        if self._run is not None:
            return
        if not self.model.placements:
            QMessageBox.information(self, "Empty cell", "Add at least one molecule.")
            return
        try:
            run = prepare_run(self.project, self.cell_x.value(),
                              self.cell_y.value(),
                              [(p.x, p.y, p.phi) for p in self.model.placements],
                              self)
        except (EngineError, OSError, ValueError) as e:
            QMessageBox.warning(self, "Cannot optimize", str(e))
            return
        self._run = run
        self._set_running(True)
        self.opt_status.setText("Optimizing…")
        run.frameReady.connect(self._apply_cell)
        run.progress.connect(self._on_opt_line)
        run.finished.connect(self._on_opt_finished)
        run.start()

    def _stop_optimize(self) -> None:
        if self._run is not None:
            self._run.stop()

    def _set_running(self, running: bool) -> None:
        self.stop_btn.setEnabled(running)
        for w in (self.optimize_btn, self.attach_btn, self.cell_x, self.cell_y,
                  self.canvas, self.table, self.add_btn, self.rem_btn,
                  self.slot_load, self.slot_clear,
                  *self.mode_buttons.values()):
            w.setEnabled(not running)

    def _apply_cell(self, cell_x: float, cell_y: float, placements: list) -> None:
        for spin, v in ((self.cell_x, cell_x), (self.cell_y, cell_y)):
            spin.blockSignals(True)
            spin.setValue(v)
            spin.blockSignals(False)
        self.canvas.cell_x, self.canvas.cell_y = cell_x, cell_y
        self.model.set_placements([Placement(round(x, 4), round(y, 4),
                                             round(phi, 2))
                                   for x, y, phi in placements])

    def _on_opt_line(self, line: str) -> None:
        if line.startswith("Density:"):
            # the density half of the line is already tracked live by the
            # info label; keep the running energy here
            self.opt_status.setText(
                f"E = {line.split('Energy:')[1].strip()} kJ/mol")
        elif line.startswith(("Cell scaling:", "Steps halved",
                              "Unit cell optimization:", "The starting cell")):
            self.opt_status.setText(line)

    def _on_opt_finished(self, ok: bool, result, message: str) -> None:
        self._run = None
        self._set_running(False)
        self._refresh_slot()
        if ok:
            self._apply_cell(*result)
            self.opt_status.setText(f"Optimized: {message}")
            self.statusMessage.emit(f"Unit cell optimized: {message}")
        else:
            self.opt_status.setText(f"Optimization {message}"
                                    if message == "stopped"
                                    else "Optimization failed")
            if message != "stopped":
                QMessageBox.warning(self, "Optimization failed", message)

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
        self.projectCellChanged.emit()
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
        self.projectCellChanged.emit()
        self.statusMessage.emit("Unit cell removed from the project")
