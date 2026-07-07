"""Tab 3 "Potentials": manage the existing pair potential of the project.
Attach a binary v2 file, inspect its header, convert an ASCII grid.

The potential file is referenced, not copied: forcefields are often
hundreds of megabytes."""

from pathlib import Path

from PySide6.QtCore import QSettings, Qt, Signal
from PySide6.QtWidgets import (QFileDialog, QGroupBox, QHBoxLayout, QLabel,
                               QMessageBox, QPushButton, QVBoxLayout, QWidget)

from ..forcefield import ForcefieldError, read_header
from ..project import Project
from .pack_dialog import PackDialog


class PotentialsTab(QWidget):
    statusMessage = Signal(str)

    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self.project = project

        root = QVBoxLayout(self)
        root.setContentsMargins(16, 16, 16, 16)
        root.setAlignment(Qt.AlignTop)

        slot = QGroupBox("Project potential")
        sl = QVBoxLayout(slot)
        self.pot_label = QLabel()
        self.pot_label.setWordWrap(True)
        self.pot_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        sl.addWidget(self.pot_label)
        row = QHBoxLayout()
        attach = QPushButton("Attach file…")
        attach.clicked.connect(self._attach)
        self.detach_btn = QPushButton("Detach")
        self.detach_btn.clicked.connect(self._detach)
        row.addWidget(attach)
        row.addWidget(self.detach_btn)
        row.addStretch(1)
        sl.addLayout(row)
        root.addWidget(slot)

        tools = QGroupBox("Convert an ASCII grid to the binary v2 format")
        tl = QVBoxLayout(tools)
        info = QLabel("Pack a text potential grid (r, θ₁, θ₂, energy) into the "
                      "compact binary format the engine reads. Optionally fold "
                      "by the molecule's rotational symmetry and subsample.")
        info.setWordWrap(True)
        info.setProperty("dim", True)
        tl.addWidget(info)
        row = QHBoxLayout()
        convert = QPushButton("Convert ASCII…")
        convert.clicked.connect(self._convert)
        row.addWidget(convert)
        row.addStretch(1)
        tl.addLayout(row)
        root.addWidget(tools)

        self.refresh()

    # -- state -------------------------------------------------------------

    def refresh(self) -> None:
        entry = self.project.potential
        if entry is None:
            self.pot_label.setText("No potential attached. Attach an existing "
                                   ".v2.bin, or convert an ASCII grid below.")
            self.detach_btn.setEnabled(False)
            return
        self.detach_btn.setEnabled(True)
        path = self.project.potential_path()
        if not path.is_file():
            self.pot_label.setText(f"<b>{entry['name']}</b><br>file missing: {path}")
            return
        try:
            header = read_header(path)
            self.pot_label.setText(f"<b>{entry['name']}</b>   ·   {path}<br>"
                                   + header.summary().replace("\n", "<br>"))
        except ForcefieldError as e:
            self.pot_label.setText(f"<b>{entry['name']}</b>   ·   {path}<br>"
                                   f"bad file: {e}")

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
