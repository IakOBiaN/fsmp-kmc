"""Tab 2 "Create potential": prepare the project model for potential
generation. The panel depends on the model kind:

- atomistic: read-only geometry, force-field optimization, partial charges
  (potential generation from an atomistic model is still planned);
- site: read-only site model overview plus a working Lennard-Jones + Coulomb
  potential generator.
"""

from pathlib import Path

from PySide6.QtCore import (QAbstractTableModel, QModelIndex, QThread, Qt,
                            Signal)
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox,
                               QDoubleSpinBox, QGroupBox, QHBoxLayout,
                               QHeaderView, QLabel, QMessageBox, QProgressBar,
                               QPushButton, QSpinBox, QSplitter, QStackedWidget,
                               QTableView, QVBoxLayout, QWidget)

from .. import chem, theme
from ..atom_table import AtomTableModel
from ..canvas import MoleculeCanvas, SiteCanvas
from ..generate import GenerationError, GridSpec, generate
from ..molecule import Molecule
from ..project import Project, safe_filename
from ..site_table import SiteTableModel
from ..sitemodel import SiteModel


def _generation_placeholder(kind: str) -> QGroupBox:
    box = QGroupBox("Potential generation")
    lay = QVBoxLayout(box)
    text = ("Sweep rigid dimer configurations with the chosen force field "
            "(MMFF94/UFF now, xtb/MOPAC later) and write a binary v2 grid.")
    label = QLabel(text + "\n\nComing next.")
    label.setWordWrap(True)
    label.setProperty("dim", True)
    lay.addWidget(label)
    badge = QLabel("planned")
    badge.setStyleSheet(f"color: {theme.GOLD}; font-weight: 600; "
                        "letter-spacing: 2px; background: transparent;")
    lay.addWidget(badge)
    return box


class AtomsChargesModel(QAbstractTableModel):
    """Read-only atoms table with a charge column."""

    COLUMNS = ["Element", "x, Å", "y, Å", "z, Å", "q, e"]

    def __init__(self, parent=None):
        super().__init__(parent)
        self._mol = Molecule()
        self._charges: list[float] | None = None

    def set_content(self, mol: Molecule, charges: list[float] | None) -> None:
        self.beginResetModel()
        self._mol = mol
        self._charges = charges
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._mol.atoms)

    def columnCount(self, parent=QModelIndex()) -> int:
        return len(self.COLUMNS)

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return None
        return self.COLUMNS[section] if orientation == Qt.Horizontal else section + 1

    def data(self, index, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return None
        a = self._mol.atoms[index.row()]
        col = index.column()
        if col == 0:
            return a.element
        if col < 4:
            return f"{(a.x, a.y, a.z)[col - 1]:.4f}"
        if self._charges is None:
            return "—"
        return f"{self._charges[index.row()]:+.4f}"


class AtomisticPage(QWidget):
    statusMessage = Signal(str)

    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self.project = project
        self._mol: Molecule | None = None
        self._preview = False
        self._charges: tuple[str, list[float]] | None = None

        splitter = QSplitter(Qt.Horizontal)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(splitter)

        self.canvas_model = AtomTableModel(self)
        self.canvas = MoleculeCanvas(self.canvas_model, read_only=True)
        self.canvas.cursorMoved.connect(
            lambda x, y: self.statusMessage.emit(f"x = {x:.2f} Å,  y = {y:.2f} Å"))
        splitter.addWidget(self.canvas)
        splitter.addWidget(self._build_side())
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 0)
        splitter.setSizes([820, 460])

    def _build_side(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(8, 8, 8, 8)

        top = QHBoxLayout()
        self.mol_label = QLabel()
        self.mol_label.setWordWrap(True)
        top.addWidget(self.mol_label, 1)
        top.addWidget(QLabel("Total charge"))
        self.charge_spin = QSpinBox()
        self.charge_spin.setRange(-8, 8)
        self.charge_spin.setValue(0)
        self.charge_spin.setToolTip("Total charge of the molecule, e")
        top.addWidget(self.charge_spin)
        layout.addLayout(top)

        if not chem.HAVE_RDKIT:
            warn = QLabel("RDKit is not installed; optimization and charges are "
                          "disabled  (pip install rdkit)")
            warn.setWordWrap(True)
            warn.setProperty("dim", True)
            layout.addWidget(warn)

        self.opt_group = QGroupBox("Geometry optimization (RDKit)")
        og = QVBoxLayout(self.opt_group)
        row = QHBoxLayout()
        row.addWidget(QLabel("Method"))
        self.opt_method = QComboBox()
        self.opt_method.addItems(chem.OPTIMIZE_METHODS)
        row.addWidget(self.opt_method)
        row.addWidget(QLabel("Max steps"))
        self.opt_steps = QSpinBox()
        self.opt_steps.setRange(10, 100000)
        self.opt_steps.setValue(500)
        row.addWidget(self.opt_steps)
        row.addStretch(1)
        self.opt_btn = QPushButton("Optimize")
        self.opt_btn.clicked.connect(self._optimize)
        row.addWidget(self.opt_btn)
        og.addLayout(row)
        self.opt_result = QLabel(" ")
        self.opt_result.setWordWrap(True)
        self.opt_result.setProperty("dim", True)
        og.addWidget(self.opt_result)
        self.apply_btn = QPushButton("Apply to project")
        self.apply_btn.setProperty("primary", True)
        self.apply_btn.setEnabled(False)
        self.apply_btn.clicked.connect(self._apply_geometry)
        og.addWidget(self.apply_btn)
        layout.addWidget(self.opt_group)

        self.chg_group = QGroupBox("Partial charges (RDKit)")
        cg = QVBoxLayout(self.chg_group)
        row = QHBoxLayout()
        row.addWidget(QLabel("Method"))
        self.chg_method = QComboBox()
        self.chg_method.addItems(chem.CHARGE_METHODS)
        row.addWidget(self.chg_method)
        row.addStretch(1)
        self.chg_btn = QPushButton("Compute")
        self.chg_btn.clicked.connect(self._compute_charges)
        row.addWidget(self.chg_btn)
        cg.addLayout(row)
        self.chg_result = QLabel(" ")
        self.chg_result.setWordWrap(True)
        self.chg_result.setProperty("dim", True)
        cg.addWidget(self.chg_result)
        self.store_btn = QPushButton("Store in project")
        self.store_btn.setEnabled(False)
        self.store_btn.clicked.connect(self._store_charges)
        cg.addWidget(self.store_btn)
        layout.addWidget(self.chg_group)

        self.atoms_model = AtomsChargesModel(self)
        table = QTableView()
        table.setModel(self.atoms_model)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        table.setSelectionBehavior(QAbstractItemView.SelectRows)
        table.setAlternatingRowColors(True)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table.verticalHeader().setDefaultSectionSize(24)
        layout.addWidget(table, 1)

        layout.addWidget(_generation_placeholder("atomistic"))
        return panel

    def refresh(self) -> None:
        self._preview = False
        self._charges = None
        self.apply_btn.setEnabled(False)
        self.store_btn.setEnabled(False)
        self.opt_result.setText(" ")
        self.chg_result.setText(" ")

        entry = self.project.model
        path = self.project.model_path(entry)
        try:
            self._mol = Molecule.load_xyz(path)
            self.mol_label.setText(f"<b>{entry['name']}</b>   ·   "
                                   f"{self._mol.formula()}")
        except (OSError, ValueError) as e:
            self._mol = None
            self.mol_label.setText(f"Cannot read {path}: {e}")

        stored = self.project.charges
        charges = None
        if (self._mol is not None and stored is not None
                and len(stored["values"]) == len(self._mol.atoms)):
            charges = stored["values"]
            self.chg_result.setText(f"stored in project: {stored['method']}")
        self.atoms_model.set_content(self._mol or Molecule(), charges)
        self.canvas_model.set_molecule(self._mol or Molecule())
        self.canvas.reset_view()

        enabled = self._mol is not None and chem.HAVE_RDKIT
        self.opt_group.setEnabled(enabled)
        self.chg_group.setEnabled(enabled)

    def showEvent(self, event):
        super().showEvent(event)
        self.canvas.reset_view()  # the viewport now has its final size

    # -- actions -----------------------------------------------------------

    def _optimize(self) -> None:
        if self._mol is None:
            return
        try:
            res = chem.optimize(self._mol, self.opt_method.currentText(),
                                self.charge_spin.value(), self.opt_steps.value())
        except chem.ChemError as e:
            QMessageBox.warning(self, "Optimization failed", str(e))
            return
        self._mol = res.molecule
        self._preview = True
        self._charges = None
        self.canvas_model.set_molecule(self._mol)
        self.atoms_model.set_content(self._mol, None)
        self.canvas.reset_view()
        state = "converged" if res.converged else "NOT converged, raise max steps"
        self.opt_result.setText(
            f"{self.opt_method.currentText()}: E {res.e_before:.2f} → "
            f"{res.e_after:.2f} kJ/mol, {state}. Not saved to the project yet.")
        self.apply_btn.setEnabled(True)
        self.store_btn.setEnabled(False)
        self.statusMessage.emit("Optimization done")

    def _apply_geometry(self) -> None:
        entry = self.project.model
        if entry is None or self._mol is None:
            return
        self.project.set_atomistic(entry["name"], self._mol)
        self._preview = False
        self.apply_btn.setEnabled(False)
        self.opt_result.setText(self.opt_result.text().replace(
            "Not saved to the project yet.", "Saved to the project."))
        self.statusMessage.emit("Optimized geometry saved to the project")

    def _compute_charges(self) -> None:
        if self._mol is None:
            return
        method = self.chg_method.currentText()
        try:
            values = chem.partial_charges(self._mol, method,
                                          self.charge_spin.value())
        except chem.ChemError as e:
            QMessageBox.warning(self, "Charge calculation failed", str(e))
            return
        self._charges = (method, values)
        self.atoms_model.set_content(self._mol, values)
        total = sum(values)
        note = (" Apply the optimized geometry before storing."
                if self._preview else "")
        self.chg_result.setText(f"{method}: sum {total:+.4f} e. "
                                f"Not stored in the project yet.{note}")
        self.store_btn.setEnabled(not self._preview)
        self.statusMessage.emit("Charges computed")

    def _store_charges(self) -> None:
        if self._charges is None or self._preview:
            return
        method, values = self._charges
        self.project.set_charges(method, values)
        self.store_btn.setEnabled(False)
        self.chg_result.setText(f"stored in project: {method}")
        self.statusMessage.emit("Charges stored in the project")


class GenerateWorker(QThread):
    progress = Signal(int, int)
    finished_ok = Signal(str)
    failed = Signal(str)

    def __init__(self, model, spec, out_path, parent=None):
        super().__init__(parent)
        self._model, self._spec, self._out = model, spec, out_path
        self._cancel = False

    def cancel(self):
        self._cancel = True

    def run(self):
        try:
            generate(self._model, self._spec, self._out,
                     progress=lambda d, t: self.progress.emit(d, t),
                     cancel=lambda: self._cancel)
        except (GenerationError, OSError, ValueError) as e:
            self.failed.emit(str(e))
            return
        if not self._cancel:
            self.finished_ok.emit(str(self._out))


class SitePage(QWidget):
    statusMessage = Signal(str)
    potentialGenerated = Signal()

    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self.project = project
        self._site_model = SiteModel()
        self._worker: GenerateWorker | None = None

        splitter = QSplitter(Qt.Horizontal)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(splitter)

        self.model = SiteTableModel(self)
        self.canvas = SiteCanvas(self.model, read_only=True)
        self.canvas.cursorMoved.connect(
            lambda x, y: self.statusMessage.emit(f"x = {x:.2f} Å,  y = {y:.2f} Å"))
        splitter.addWidget(self.canvas)
        splitter.addWidget(self._build_side())
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 0)
        splitter.setSizes([760, 520])

    def _build_side(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(8, 8, 8, 8)

        self.title = QLabel()
        self.title.setWordWrap(True)
        layout.addWidget(self.title)

        table = QTableView()
        table.setModel(self.model)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        table.setSelectionBehavior(QAbstractItemView.SelectRows)
        table.setAlternatingRowColors(True)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table.verticalHeader().setDefaultSectionSize(24)
        layout.addWidget(table, 1)

        layout.addWidget(self._build_generation())
        return panel

    def _build_generation(self) -> QGroupBox:
        box = QGroupBox("Potential generation  (Lennard-Jones + Coulomb)")
        v = QVBoxLayout(box)

        def spin(lo, hi, val, step, dec, suffix):
            s = QDoubleSpinBox()
            s.setRange(lo, hi)
            s.setValue(val)
            s.setSingleStep(step)
            s.setDecimals(dec)
            s.setSuffix(suffix)
            return s

        self.r_min = spin(0.1, 100, 7.6, 0.1, 3, " Å")
        self.r_max = spin(1, 200, 30.0, 1, 2, " Å")
        self.dr = spin(0.005, 5, 0.05, 0.01, 3, " Å")
        self.da = spin(0.5, 60, 5.0, 0.5, 1, "°")
        self.fold = QSpinBox()
        self.fold.setRange(0, 360)
        self.fold.setValue(0)
        self.fold.setSpecialValueText("none")
        self.fold.setSuffix("°")
        self.fold.setToolTip("Rotational symmetry of the molecule (120 for TMA); "
                             "0 = no folding. Must divide 360.")
        self.use_float = QCheckBox("float (half the size)")
        self.use_float.setChecked(True)

        r1 = QHBoxLayout()
        for lab, w in (("r from", self.r_min), ("to", self.r_max), ("step", self.dr)):
            r1.addWidget(QLabel(lab))
            r1.addWidget(w)
        r1.addStretch(1)
        v.addLayout(r1)
        r2 = QHBoxLayout()
        r2.addWidget(QLabel("angle step"))
        r2.addWidget(self.da)
        r2.addWidget(QLabel("fold"))
        r2.addWidget(self.fold)
        r2.addWidget(self.use_float)
        r2.addStretch(1)
        v.addLayout(r2)

        row = QHBoxLayout()
        self.gen_btn = QPushButton("Generate…")
        self.gen_btn.setProperty("primary", True)
        self.gen_btn.clicked.connect(self._generate)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setEnabled(False)
        self.cancel_btn.clicked.connect(self._cancel)
        row.addWidget(self.gen_btn)
        row.addWidget(self.cancel_btn)
        row.addStretch(1)
        v.addLayout(row)

        self.progress = QProgressBar()
        self.progress.setVisible(False)
        v.addWidget(self.progress)
        self.gen_status = QLabel(" ")
        self.gen_status.setWordWrap(True)
        self.gen_status.setProperty("dim", True)
        v.addWidget(self.gen_status)
        return box

    def refresh(self) -> None:
        entry = self.project.model
        try:
            sm = SiteModel.load(self.project.model_path(entry))
        except (OSError, ValueError) as e:
            self._site_model = SiteModel()
            self.model.set_site_model(SiteModel())
            self.title.setText(f"Cannot read site model: {e}")
            return
        self._site_model = sm
        self.model.set_site_model(sm)
        self.title.setText(f"<b>{entry['name']}</b>   ·   {sm.summary()}   ·   "
                           f"total charge {sm.total_charge():+.4f} e")
        # prefill r_min just outside the hard core
        r0 = max((s.r0 for s in sm.sites if s.is_lj), default=0.0)
        if r0 > 0:
            self.r_min.setValue(round(r0 + 0.02, 2))
        self.canvas.reset_view()

    def showEvent(self, event):
        super().showEvent(event)
        self.canvas.reset_view()

    # -- generation --------------------------------------------------------

    def _generate(self) -> None:
        if self._worker is not None:
            return
        if not self._site_model.sites:
            QMessageBox.information(self, "No model", "The site model is empty.")
            return
        spec = GridSpec(self.r_min.value(), self.r_max.value(), self.dr.value(),
                        self.da.value(), self.fold.value() or 360.0,
                        self.use_float.isChecked())
        if spec.r_max <= spec.r_min or spec.n_dist < 2:
            QMessageBox.warning(self, "Bad grid", "Check the r range and step.")
            return
        entry = self.project.model
        out = self.project.root / f"{safe_filename(entry['name'])}.v2.bin"
        if out.exists():
            if QMessageBox.question(self, "Overwrite", f"{out.name} exists. "
                                    "Overwrite?") != QMessageBox.Yes:
                return
        self._out = out
        self.progress.setVisible(True)
        self.progress.setValue(0)
        self.gen_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.gen_status.setText(f"Generating {spec.n_dist} distances…")
        self._worker = GenerateWorker(self._site_model, spec, out, self)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished_ok.connect(self._on_done)
        self._worker.failed.connect(self._on_failed)
        self._worker.start()

    def _on_progress(self, done: int, total: int) -> None:
        self.progress.setMaximum(total)
        self.progress.setValue(done)

    def _cancel(self) -> None:
        if self._worker is not None:
            self._worker.cancel()
            self.gen_status.setText("Cancelling…")

    def _cleanup_worker(self) -> None:
        self._worker = None
        self.progress.setVisible(False)
        self.gen_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)

    def _on_done(self, path: str) -> None:
        self._cleanup_worker()
        self.gen_status.setText(f"Done: {Path(path).name}")
        self.statusMessage.emit(f"Potential written to {path}")
        if QMessageBox.question(self, "Attach potential",
                                "Attach the generated potential to the project?"
                                ) == QMessageBox.Yes:
            name = Path(path).name.removesuffix(".v2.bin").removesuffix(".bin")
            self.project.set_potential(name, path)
            self.potentialGenerated.emit()

    def _on_failed(self, msg: str) -> None:
        self._cleanup_worker()
        self.gen_status.setText(f"Failed: {msg}")
        QMessageBox.warning(self, "Generation failed", msg)


class CreatePotentialTab(QWidget):
    statusMessage = Signal(str)

    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self.project = project

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)

        self.stack = QStackedWidget()
        self.empty_page = QLabel("No project model yet. Add an atomistic "
                                 "molecule or a site model on the Molecule "
                                 "model tab first.")
        self.empty_page.setAlignment(Qt.AlignCenter)
        self.empty_page.setProperty("dim", True)
        self.atomistic_page = AtomisticPage(project)
        self.site_page = SitePage(project)
        for w in (self.empty_page, self.atomistic_page, self.site_page):
            self.stack.addWidget(w)
        self.atomistic_page.statusMessage.connect(self.statusMessage)
        self.site_page.statusMessage.connect(self.statusMessage)
        layout.addWidget(self.stack)

        self.refresh()

    def refresh(self) -> None:
        kind = self.project.model_kind
        if kind == "atomistic":
            self.atomistic_page.refresh()
            self.stack.setCurrentWidget(self.atomistic_page)
        elif kind == "site":
            self.site_page.refresh()
            self.stack.setCurrentWidget(self.site_page)
        else:
            self.stack.setCurrentWidget(self.empty_page)
