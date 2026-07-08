"""Tab 2 "Create potential": prepare the project model for potential
generation. The panel depends on the model kind:

- atomistic: read-only geometry; generation via ANI-2x (torchani) is planned,
  so the panel is a placeholder for now;
- site: read-only site model overview plus a working Lennard-Jones + Coulomb
  potential generator.
"""

from pathlib import Path

from PySide6.QtCore import QThread, Qt, Signal
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QDoubleSpinBox,
                               QGroupBox, QHBoxLayout, QHeaderView, QLabel,
                               QMessageBox, QProgressBar, QPushButton, QSpinBox,
                               QSplitter, QStackedWidget, QTableView,
                               QVBoxLayout, QWidget)

from .. import theme
from ..atom_table import AtomTableModel
from ..canvas import MoleculeCanvas, SiteCanvas
from ..generate import GenerationError, GridSpec, generate
from ..molecule import Molecule
from ..project import Project, safe_filename
from ..site_table import SiteTableModel
from ..sitemodel import SiteModel


class AtomisticPage(QWidget):
    """Read-only view of the atomistic project molecule. Potential generation
    from an atomistic model is planned via ANI-2x (torchani) and not built
    yet, so the panel is a placeholder."""

    statusMessage = Signal(str)

    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self.project = project

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
        splitter.setSizes([820, 440])

    def _build_side(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(8, 8, 8, 8)

        self.title = QLabel()
        self.title.setWordWrap(True)
        layout.addWidget(self.title)

        table = QTableView()
        table.setModel(self.canvas_model)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        table.setSelectionBehavior(QAbstractItemView.SelectRows)
        table.setAlternatingRowColors(True)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table.verticalHeader().setDefaultSectionSize(24)
        layout.addWidget(table, 1)

        box = QGroupBox("Potential generation")
        v = QVBoxLayout(box)
        info = QLabel("Potential generation from an atomistic molecule will use "
                      "the ANI-2x machine-learning potential (torchani).")
        info.setWordWrap(True)
        info.setProperty("dim", True)
        v.addWidget(info)
        badge = QLabel("planned")
        badge.setStyleSheet(f"color: {theme.GOLD}; font-weight: 600; "
                            "letter-spacing: 2px; background: transparent;")
        v.addWidget(badge)
        layout.addWidget(box)
        return panel

    def refresh(self) -> None:
        entry = self.project.atomistic
        try:
            mol = Molecule.load_xyz(self.project.model_path(entry))
        except (OSError, ValueError) as e:
            self.canvas_model.set_molecule(Molecule())
            self.title.setText(f"Cannot read molecule: {e}")
            return
        self.canvas_model.set_molecule(mol)
        self.title.setText(f"<b>{entry['name']}</b>   ·   {mol.formula()}   ·   "
                           f"{len(mol.atoms)} atoms")
        self.canvas.reset_view()

    def showEvent(self, event):
        super().showEvent(event)
        self.canvas.reset_view()


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
        entry = self.project.site
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
        entry = self.project.site
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
        # the site model drives generation; a lone atomistic model falls back
        # to the (planned) ANI-2x page
        kind = self.project.generation_kind
        if kind == "site":
            self.site_page.refresh()
            self.stack.setCurrentWidget(self.site_page)
        elif kind == "atomistic":
            self.atomistic_page.refresh()
            self.stack.setCurrentWidget(self.atomistic_page)
        else:
            self.stack.setCurrentWidget(self.empty_page)
