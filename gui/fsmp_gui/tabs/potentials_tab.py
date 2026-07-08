"""Tab 3 "Potentials": manage the project's pair potential (attach a binary
v2 file, inspect its header, convert an ASCII grid) and browse its values with
a two-molecule viewer.

The potential file is referenced, not copied: forcefields are often hundreds
of megabytes, so the viewer reads values on demand from a memory map."""

from pathlib import Path

from PySide6.QtCore import QSettings, Qt, Signal
from PySide6.QtWidgets import (QDoubleSpinBox, QFileDialog, QGroupBox,
                               QHBoxLayout, QLabel, QMessageBox, QPushButton,
                               QVBoxLayout, QWidget)

from ..canvas import site_color
from ..elements import covalent_radius, element_color
from ..forcefield import ForcefieldError, ForcefieldGrid, read_header
from ..molecule import Molecule
from ..project import Project
from ..sitemodel import SiteModel
from ..potential_viewer import TwoMoleculeView
from .pack_dialog import PackDialog


class PotentialsTab(QWidget):
    statusMessage = Signal(str)

    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self.project = project
        self._grid: ForcefieldGrid | None = None

        root = QVBoxLayout(self)
        root.setContentsMargins(16, 16, 16, 16)

        slot = QGroupBox("Project potential")
        sl = QVBoxLayout(slot)
        self.pot_label = QLabel()
        self.pot_label.setWordWrap(True)
        self.pot_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        sl.addWidget(self.pot_label)
        row = QHBoxLayout()
        attach = QPushButton("Attach file…")
        attach.clicked.connect(self._attach)
        convert = QPushButton("Convert ASCII…")
        convert.clicked.connect(self._convert)
        self.detach_btn = QPushButton("Detach")
        self.detach_btn.clicked.connect(self._detach)
        row.addWidget(attach)
        row.addWidget(convert)
        row.addWidget(self.detach_btn)
        row.addStretch(1)
        sl.addLayout(row)
        root.addWidget(slot)

        root.addWidget(self._build_viewer(), 1)
        self.refresh()

    def _build_viewer(self) -> QGroupBox:
        box = QGroupBox("Potential viewer")
        v = QVBoxLayout(box)
        self.view = TwoMoleculeView()
        v.addWidget(self.view, 1)

        controls = QHBoxLayout()
        self.r_spin = QDoubleSpinBox()
        self.r_spin.setDecimals(2)
        self.r_spin.setSingleStep(0.1)
        self.r_spin.setSuffix(" Å")
        self.t1_spin = QDoubleSpinBox()
        self.t2_spin = QDoubleSpinBox()
        for s in (self.t1_spin, self.t2_spin):
            s.setRange(0.0, 360.0)
            s.setSingleStep(5.0)
            s.setDecimals(1)
            s.setSuffix("°")
            s.setWrapping(True)
        for lab, s in (("distance r", self.r_spin), ("angle θ₁ (A)", self.t1_spin),
                       ("angle θ₂ (B)", self.t2_spin)):
            controls.addWidget(QLabel(lab))
            controls.addWidget(s)
        controls.addStretch(1)
        self.energy_label = QLabel()
        self.energy_label.setProperty("gold", True)
        self.energy_label.setStyleSheet("font-size: 12pt; font-weight: 600; "
                                        "background: transparent;")
        controls.addWidget(self.energy_label)
        v.addLayout(controls)

        for s in (self.r_spin, self.t1_spin, self.t2_spin):
            s.valueChanged.connect(self._update_view)
        return box

    # -- state -------------------------------------------------------------

    def refresh(self) -> None:
        self._grid = None
        entry = self.project.potential
        if entry is None:
            self.pot_label.setText("No potential attached. Attach an existing "
                                   ".v2.bin, or convert an ASCII grid.")
            self.detach_btn.setEnabled(False)
        else:
            self.detach_btn.setEnabled(True)
            path = self.project.potential_path()
            if not path.is_file():
                self.pot_label.setText(f"<b>{entry['name']}</b><br>file missing: {path}")
            else:
                try:
                    header = read_header(path)
                    self.pot_label.setText(f"<b>{entry['name']}</b>   ·   {path}<br>"
                                           + header.summary().replace("\n", "<br>"))
                    self._grid = ForcefieldGrid(header)
                except (ForcefieldError, OSError, ValueError) as e:
                    self.pot_label.setText(f"<b>{entry['name']}</b>   ·   {path}<br>"
                                           f"bad file: {e}")

        self.view.set_glyph(self._build_glyph())
        self._configure_spins()
        self._update_view()

    def _configure_spins(self) -> None:
        enabled = self._grid is not None
        for s in (self.r_spin, self.t1_spin, self.t2_spin):
            s.blockSignals(True)
            s.setEnabled(enabled)
        if enabled:
            info = self._grid.info
            self.r_spin.setRange(round(info.min_dist, 2), round(info.r_max, 2))
            self.r_spin.setValue(round(min(info.min_dist + 4.0, info.r_max), 2))
            self.t1_spin.setValue(0.0)
            self.t2_spin.setValue(60.0)
        for s in (self.r_spin, self.t1_spin, self.t2_spin):
            s.blockSignals(False)

    def _build_glyph(self):
        entry = self.project.model
        if entry is None:
            return None
        try:
            path = self.project.model_path(entry)
            if entry["kind"] == "atomistic":
                mol = Molecule.load_xyz(path)
                return [(a.x, a.y, element_color(a.element),
                         min(max(0.4 * covalent_radius(a.element) + 0.15, 0.3), 0.8))
                        for a in mol.atoms]
            sm = SiteModel.load(path)
            return [(s.x, s.y, site_color(s.q).name(),
                     0.32 if s.sigma <= 0 else min(max(0.2 * s.sigma + 0.1, 0.32), 0.9))
                    for s in sm.sites]
        except (OSError, ValueError):
            return None

    def _update_view(self) -> None:
        r = self.r_spin.value()
        t1, t2 = self.t1_spin.value(), self.t2_spin.value()
        self.view.set_config(r, t1, t2)
        if self._grid is None:
            self.energy_label.setText("attach a potential to read energies")
            return
        e = self._grid.energy_at(r, t1, t2)
        if e is None:
            self.energy_label.setText("hard core (r below the minimum)")
        else:
            self.energy_label.setText(f"U = {e:,.1f} J/mol   ({e / 1000:.3f} kJ/mol)")

    def showEvent(self, event):
        super().showEvent(event)
        self._update_view()

    # -- actions -----------------------------------------------------------

    def _attach(self) -> None:
        start = QSettings().value("potential/last_dir", "")
        path, _ = QFileDialog.getOpenFileName(
            self, "Attach potential", start,
            "FSMP forcefields (*.bin);;All files (*)")
        if not path:
            return
        QSettings().setValue("potential/last_dir", str(Path(path).parent))
        try:
            read_header(path)
        except ForcefieldError as e:
            QMessageBox.warning(self, "Not a v2 forcefield", f"{path}\n\n{e}")
            return
        name = Path(path).name.removesuffix(".v2.bin").removesuffix(".bin")
        self.project.set_potential(name, path)
        self.refresh()
        self.statusMessage.emit(f"Potential '{name}' attached to the project")

    def _convert(self) -> None:
        dialog = PackDialog(self)
        dialog.exec()
        if dialog.success and dialog.output_path:
            answer = QMessageBox.question(
                self, "Attach converted potential",
                "Attach the converted potential to the project?")
            if answer == QMessageBox.Yes:
                name = Path(dialog.output_path).name
                name = name.removesuffix(".v2.bin").removesuffix(".bin")
                self.project.set_potential(name, dialog.output_path)
                self.refresh()
                self.statusMessage.emit(f"Potential '{name}' attached to the project")

    def _detach(self) -> None:
        self.project.clear_potential()
        self.refresh()
        self.statusMessage.emit("Potential detached (the file is kept)")
