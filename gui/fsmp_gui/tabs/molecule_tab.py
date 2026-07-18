"""Atomistic molecule editor (subtab of "Molecule model").

Build, edit, load and save an atomistic geometry; "Use in project" makes it
the project's atomistic model (used for visualization).
"""

from pathlib import Path

from PySide6.QtCore import QSettings, Qt, Signal
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QButtonGroup,
                               QComboBox, QDoubleSpinBox, QFileDialog, QFrame,
                               QHBoxLayout, QHeaderView, QInputDialog, QLabel,
                               QMessageBox, QPushButton, QSplitter, QTableView,
                               QToolButton, QVBoxLayout, QWidget)

from .. import mmff
from ..atom_table import AtomTableModel
from ..canvas import Mode, MoleculeCanvas
from ..elements import COMMON, normalize_symbol
from ..molecule import (Atom, Molecule, aimed_at_x, centroid,
                        connected_components, rotated, translated)
from ..project import Project


class MoleculeTab(QWidget):
    statusMessage = Signal(str)
    projectModelChanged = Signal()

    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self.project = project
        self.model = AtomTableModel(self)
        self.model.changed.connect(self._on_model_changed)
        self._dirty = False
        self._current_name: str | None = None

        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.addLayout(self._build_toolbar())

        splitter = QSplitter(Qt.Horizontal)
        root.addWidget(splitter, 1)

        self.canvas = MoleculeCanvas(self.model)
        self.canvas.cursorMoved.connect(
            lambda x, y: self.statusMessage.emit(f"x = {x:.2f} Å,  y = {y:.2f} Å"))
        self.canvas.addRequested.connect(self._add_atom_from_canvas)
        splitter.addWidget(self.canvas)
        splitter.addWidget(self._build_side_panel())
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 0)
        splitter.setSizes([900, 380])

        self._set_mode(Mode.SELECT)
        self.refresh_project_model()
        # the attached model goes straight into the editor
        entry = project.atomistic
        if entry is not None:
            self._load_file(str(project.model_path(entry)), name=entry["name"])
        self._update_info()
        self._fit_on_show = True

    def showEvent(self, event):
        super().showEvent(event)
        if self._fit_on_show:   # the pre-show fit had no real viewport size
            self.canvas.reset_view()
            self._fit_on_show = False

    # -- construction ------------------------------------------------------

    def _build_toolbar(self) -> QHBoxLayout:
        bar = QHBoxLayout()
        bar.setSpacing(6)

        for text, slot in (("New", self.new_molecule),
                           ("Open…", self.open_molecule),
                           ("Save as…", self.save_molecule)):
            b = QPushButton(text)
            b.clicked.connect(slot)
            bar.addWidget(b)

        self.optimize_btn = QPushButton("Optimize (MMFF94)")
        self.optimize_btn.clicked.connect(self.optimize_geometry)
        if mmff.rdkit_available():
            self.optimize_btn.setToolTip("Relax the geometry with the MMFF94 "
                                         "force field (kept planar)")
        else:
            self.optimize_btn.setEnabled(False)
            self.optimize_btn.setToolTip("Install RDKit (pip install rdkit) to "
                                         "optimize the geometry")
        bar.addWidget(self.optimize_btn)

        bar.addSpacing(16)

        center = QPushButton("Center")
        center.setToolTip("Move the molecule centre to the origin, the point "
                          "the engine rotates the molecule about")
        center.clicked.connect(self.center_molecule)
        bar.addWidget(center)

        ccw = QToolButton()
        ccw.setText("↺")
        ccw.setToolTip("Rotate the whole molecule counterclockwise by the step")
        ccw.clicked.connect(lambda: self.rotate_molecule(+1.0))
        bar.addWidget(ccw)
        self.rotate_step = QDoubleSpinBox()
        self.rotate_step.setRange(0.1, 180.0)
        self.rotate_step.setValue(15.0)
        self.rotate_step.setDecimals(1)
        self.rotate_step.setSingleStep(5.0)
        self.rotate_step.setSuffix("°")
        self.rotate_step.setToolTip("Rotation step")
        bar.addWidget(self.rotate_step)
        cw = QToolButton()
        cw.setText("↻")
        cw.setToolTip("Rotate the whole molecule clockwise by the step")
        cw.clicked.connect(lambda: self.rotate_molecule(-1.0))
        bar.addWidget(cw)

        aim = QPushButton("Point at +x")
        aim.setToolTip("Rotate the molecule rigidly so the selected atom lies "
                       "on the positive x axis, the zero-angle direction")
        aim.clicked.connect(self.aim_selected_at_x)
        bar.addWidget(aim)

        bar.addSpacing(16)

        self.mode_group = QButtonGroup(self)
        self.mode_buttons = {}
        for mode, text, tip in ((Mode.SELECT, "Select", "Select and drag atoms"),
                                (Mode.ADD, "Add atom", "Click on empty space to add an atom"),
                                (Mode.DELETE, "Delete", "Click an atom to delete it")):
            b = QToolButton()
            b.setText(text)
            b.setToolTip(tip)
            b.setCheckable(True)
            b.clicked.connect(lambda _=False, m=mode: self._set_mode(m))
            self.mode_group.addButton(b)
            self.mode_buttons[mode] = b
            bar.addWidget(b)

        self.element_box = QComboBox()
        self.element_box.setEditable(True)
        self.element_box.addItems(COMMON)
        self.element_box.setToolTip("Element for new atoms (type any symbol)")
        self.element_box.setFixedWidth(72)
        bar.addWidget(self.element_box)

        bar.addSpacing(16)
        reset = QPushButton("Fit view")
        reset.clicked.connect(self.canvas_reset)
        bar.addWidget(reset)

        bar.addStretch(1)

        attach = QPushButton("Use in project")
        attach.setProperty("primary", True)
        attach.setToolTip("Save this molecule as the project molecule")
        attach.clicked.connect(self.use_in_project)
        bar.addWidget(attach)
        return bar

    def _build_side_panel(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(8, 0, 0, 0)

        caption = QLabel("Atoms")
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
        add = QPushButton("Add row")
        add.clicked.connect(self._add_row)
        rem = QPushButton("Remove selected")
        rem.clicked.connect(self._remove_selected_rows)
        row.addWidget(add)
        row.addWidget(rem)
        row.addStretch(1)
        layout.addLayout(row)

        self.info = QLabel()
        self.info.setProperty("dim", True)
        layout.addWidget(self.info)

        layout.addSpacing(10)
        caption2 = QLabel("Atomistic model in project  (visualization)")
        caption2.setProperty("dim", True)
        layout.addWidget(caption2)

        frame = QFrame()
        frame.setProperty("panel", True)
        frame_layout = QVBoxLayout(frame)
        self.slot_label = QLabel()
        self.slot_label.setWordWrap(True)
        frame_layout.addWidget(self.slot_label)
        slot_row = QHBoxLayout()
        self.slot_open = QPushButton("Open in editor")
        self.slot_open.clicked.connect(self._open_project_molecule)
        self.slot_clear = QPushButton("Remove")
        self.slot_clear.clicked.connect(self._clear_project_molecule)
        slot_row.addWidget(self.slot_open)
        slot_row.addWidget(self.slot_clear)
        slot_row.addStretch(1)
        frame_layout.addLayout(slot_row)
        layout.addWidget(frame)
        return panel

    # -- state helpers -----------------------------------------------------

    def _on_model_changed(self) -> None:
        self._dirty = True
        self._update_info()

    def _update_info(self) -> None:
        mol = self.model.molecule
        formula = mol.formula() or "empty"
        self.info.setText(f"{formula}   ·   {len(mol.atoms)} atoms")

    def _confirm_discard(self) -> bool:
        if not self._dirty or not self.model.molecule.atoms:
            return True
        answer = QMessageBox.question(
            self, "Unsaved changes",
            "The current molecule has unsaved changes. Discard them?")
        return answer == QMessageBox.Yes

    def _set_mode(self, mode: Mode) -> None:
        self.mode_buttons[mode].setChecked(True)
        self.canvas.set_mode(mode)

    def _picked_element(self) -> str:
        symbol = normalize_symbol(self.element_box.currentText())
        return symbol if symbol else "C"

    def canvas_reset(self) -> None:
        self.canvas.reset_view()

    # -- editing actions ---------------------------------------------------

    def _add_atom_from_canvas(self, x: float, y: float) -> None:
        self.canvas.add_atom_at(self._picked_element(), x, y)

    def _add_row(self) -> None:
        self.model.add_atom(Atom(self._picked_element(), 0.0, 0.0))

    def _remove_selected_rows(self) -> None:
        rows = [i.row() for i in self.table.selectionModel().selectedRows()]
        self.model.remove_rows(rows)

    def _connected_or_warn(self, mol: Molecule) -> bool:
        """A molecule model must be whole. Warn and return False when it has
        atoms or parts that are not bonded to the rest."""
        if len(mol.atoms) < 2:
            return True
        parts = connected_components(mol.atoms)
        if len(parts) == 1:
            return True
        parts.sort(key=len, reverse=True)
        stray = sorted(i for part in parts[1:] for i in part)
        names = ", ".join(f"{mol.atoms[i].element}{i + 1}" for i in stray[:10])
        if len(stray) > 10:
            names += f", … (+{len(stray) - 10})"
        QMessageBox.warning(
            self, "Molecule is not connected",
            f"The molecule has {len(parts)} separate parts: some atoms are not "
            f"bonded to the main structure ({names}).\n\n"
            "A molecule must be a single connected structure. Move the stray "
            "atoms into bonding range or remove them.")
        return False

    # -- whole-molecule geometry ---------------------------------------------

    def _apply_geometry(self, atoms: list[Atom]) -> None:
        self.model.set_molecule(Molecule(atoms, self.model.molecule.comment))

    def center_molecule(self) -> None:
        """Shift the molecule so its centre (the engine's rotation point)
        sits at the origin."""
        atoms = self.model.molecule.atoms
        if not atoms:
            return
        cx, cy = centroid(atoms)
        self._apply_geometry(translated(atoms, -cx, -cy))
        self.canvas.reset_view()
        self.statusMessage.emit("Molecule centred on its centroid")

    def rotate_molecule(self, direction: float) -> None:
        """Rotate the whole molecule about the origin by the toolbar step."""
        atoms = self.model.molecule.atoms
        if not atoms:
            return
        deg = direction * self.rotate_step.value()
        self._apply_geometry(rotated(atoms, deg))
        self.statusMessage.emit(f"Rotated by {deg:+.1f}°")

    def aim_selected_at_x(self) -> None:
        """Turn the molecule rigidly so the selected atom points along +x,
        making the zero orientation angle a chosen feature of the molecule."""
        rows = self.canvas.selected_rows()
        if not rows:
            rows = [i.row() for i in self.table.selectionModel().selectedRows()]
        if len(rows) != 1:
            QMessageBox.information(
                self, "Point at +x",
                "Select exactly one atom, on the canvas or in the table; the "
                "molecule will turn so that atom points along +x.")
            return
        atoms = self.model.molecule.atoms
        try:
            aimed = aimed_at_x(atoms, rows[0])
        except ValueError as e:
            QMessageBox.information(self, "Point at +x", str(e))
            return
        label = f"{atoms[rows[0]].element}{rows[0] + 1}"
        self._apply_geometry(aimed)
        self.statusMessage.emit(f"{label} now points along +x")

    def optimize_geometry(self) -> None:
        """Relax the current geometry with MMFF94 (UFF fallback). The result is
        planar and recentred on the molecule centre, ready for the engine."""
        mol = self.model.molecule
        if not mol.atoms:
            QMessageBox.information(self, "Nothing to optimize",
                                    "Add some atoms first.")
            return
        if not self._connected_or_warn(mol):
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            relaxed, report = mmff.optimize_molecule(mol)
        except mmff.MMFFError as e:
            QMessageBox.warning(self, "Optimization failed", str(e))
            return
        finally:
            QApplication.restoreOverrideCursor()
        self.model.set_molecule(relaxed)
        self._dirty = True
        self.canvas.reset_view()
        self._update_info()
        self.statusMessage.emit(f"Optimized  ·  {report.summary()}")

    # -- file actions ------------------------------------------------------

    def _last_dir(self) -> str:
        return QSettings().value("molecule/last_dir", str(Path.home()))

    def _remember_dir(self, path: str) -> None:
        QSettings().setValue("molecule/last_dir", str(Path(path).parent))

    def new_molecule(self) -> None:
        if not self._confirm_discard():
            return
        self.model.set_molecule(Molecule())
        self._current_name = None
        self._dirty = False
        self.statusMessage.emit("New molecule")

    def open_molecule(self) -> None:
        if not self._confirm_discard():
            return
        path, _ = QFileDialog.getOpenFileName(
            self, "Open molecule", self._last_dir(), "XYZ files (*.xyz);;All files (*)")
        if not path:
            return
        self._load_file(path)

    def _load_file(self, path: str, name: str | None = None) -> None:
        try:
            mol = Molecule.load_xyz(path)
        except (OSError, ValueError) as e:
            QMessageBox.warning(self, "Cannot open molecule", f"{path}\n\n{e}")
            return
        self.model.set_molecule(mol)
        self._current_name = name if name is not None else Path(path).stem
        self._dirty = False
        self._remember_dir(path)
        self.canvas.reset_view()
        self.statusMessage.emit(f"Opened {path}")

    def save_molecule(self) -> None:
        if not self.model.molecule.atoms:
            QMessageBox.information(self, "Nothing to save", "The molecule is empty.")
            return
        suggestion = str(Path(self._last_dir())
                         / f"{self._current_name or self.model.molecule.formula()}.xyz")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save molecule", suggestion, "XYZ files (*.xyz)")
        if not path:
            return
        try:
            self.model.molecule.save_xyz(path)
        except OSError as e:
            QMessageBox.warning(self, "Cannot save molecule", str(e))
            return
        self._current_name = Path(path).stem
        self._dirty = False
        self._remember_dir(path)
        self.statusMessage.emit(f"Saved {path}")

    # -- project integration -----------------------------------------------

    def refresh_project_model(self) -> None:
        entry = self.project.atomistic
        if entry is None:
            self.slot_label.setText("Not set yet. Build a molecule and press "
                                    "'Use in project'. Used for visualization.")
            self.slot_open.setEnabled(False)
            self.slot_clear.setEnabled(False)
        else:
            self.slot_label.setText(f"{entry['name']}   ·   {entry['file']}")
            self.slot_open.setEnabled(True)
            self.slot_clear.setEnabled(True)

    def use_in_project(self) -> None:
        mol = self.model.molecule
        if not mol.atoms:
            QMessageBox.information(self, "Empty molecule",
                                    "Add some atoms before using the molecule "
                                    "in the project.")
            return
        if not self._connected_or_warn(mol):
            return
        default = self._current_name or mol.formula()
        name, ok = QInputDialog.getText(self, "Use molecule in project",
                                        "Molecule name:", text=default)
        if not ok or not name.strip():
            return
        name = name.strip()
        old = self.project.atomistic
        if old is not None and QMessageBox.question(
                self, "Replace atomistic model",
                f"The atomistic model '{old['name']}' will be replaced. Continue?"
                ) != QMessageBox.Yes:
            return
        self.project.set_atomistic(name, mol)
        self._current_name = name
        self._dirty = False
        self.refresh_project_model()
        self.projectModelChanged.emit()
        self.statusMessage.emit(f"'{name}' is now the atomistic model")

    def _open_project_molecule(self) -> None:
        entry = self.project.atomistic
        if entry is None:
            return
        if not self._confirm_discard():
            return
        self._load_file(str(self.project.model_path(entry)), name=entry["name"])

    def _clear_project_molecule(self) -> None:
        entry = self.project.atomistic
        if entry is None:
            return
        answer = QMessageBox.question(
            self, "Remove atomistic model",
            f"Remove '{entry['name']}' from the project?\n"
            "The file inside the project folder will be deleted.")
        if answer != QMessageBox.Yes:
            return
        self.project.clear_atomistic()
        self.refresh_project_model()
        self.projectModelChanged.emit()
        self.statusMessage.emit("Atomistic model removed")
