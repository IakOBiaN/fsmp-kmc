"""Dialog that runs the pack_forcefield converter (ASCII grid -> binary v2)
as a subprocess and streams its output.

CLI of the tool (tools/pack_forcefield.cpp):
    pack [--float] [--stride N] <input.dat> <output.bin> [fold_deg]
"""

from pathlib import Path

from PySide6.QtCore import QProcess, QSettings
from PySide6.QtWidgets import (QCheckBox, QDialog, QFileDialog, QGridLayout,
                               QHBoxLayout, QLabel, QLineEdit, QMessageBox,
                               QPlainTextEdit, QPushButton, QSpinBox,
                               QVBoxLayout)


def find_default_pack() -> str:
    saved = QSettings().value("pack/exe", "")
    if saved and Path(saved).is_file():
        return saved
    repo = Path(__file__).resolve().parents[3]
    for name in ("pack.exe", "pack.out", "tools/pack.exe"):
        candidate = repo / name
        if candidate.is_file():
            return str(candidate)
    return ""


class PackDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Convert ASCII grid to binary v2")
        self.setMinimumSize(680, 520)
        self.process: QProcess | None = None
        self.success = False
        self.output_path = ""

        grid = QGridLayout()
        grid.setColumnStretch(1, 1)

        grid.addWidget(QLabel("ASCII grid"), 0, 0)
        self.input_edit = QLineEdit()
        grid.addWidget(self.input_edit, 0, 1)
        b = QPushButton("Browse…")
        b.clicked.connect(self._browse_input)
        grid.addWidget(b, 0, 2)

        grid.addWidget(QLabel("Output .bin"), 1, 0)
        self.output_edit = QLineEdit()
        grid.addWidget(self.output_edit, 1, 1)
        b = QPushButton("Browse…")
        b.clicked.connect(self._browse_output)
        grid.addWidget(b, 1, 2)

        options = QHBoxLayout()
        self.float_box = QCheckBox("float (half the size, same physics)")
        self.float_box.setChecked(True)
        options.addWidget(self.float_box)
        options.addSpacing(16)
        options.addWidget(QLabel("Fold, °"))
        self.fold_spin = QSpinBox()
        self.fold_spin.setRange(0, 180)
        self.fold_spin.setSpecialValueText("none")
        self.fold_spin.setToolTip("Rotational symmetry of the molecule "
                                  "(120 for C3, 180 for C2); must divide 360")
        options.addWidget(self.fold_spin)
        options.addSpacing(16)
        options.addWidget(QLabel("Step"))
        self.stride_spin = QSpinBox()
        self.stride_spin.setRange(1, 16)
        self.stride_spin.setToolTip("Keep every Nth grid point, coarsening the "
                                    "grid: step 2 doubles dr and dθ (1 = full grid)")
        options.addWidget(self.stride_spin)
        options.addStretch(1)
        grid.addLayout(options, 2, 0, 1, 3)

        grid.addWidget(QLabel("pack executable"), 3, 0)
        self.exe_edit = QLineEdit(find_default_pack())
        self.exe_edit.setPlaceholderText(
            "pack.exe from the release zip, or build with 'make pack'")
        grid.addWidget(self.exe_edit, 3, 1)
        b = QPushButton("Browse…")
        b.clicked.connect(self._browse_exe)
        grid.addWidget(b, 3, 2)

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setPlaceholderText("Converter output will appear here. Large "
                                    "grids take many minutes to convert.")

        buttons = QHBoxLayout()
        buttons.addStretch(1)
        self.run_btn = QPushButton("Run")
        self.run_btn.setProperty("primary", True)
        self.run_btn.clicked.connect(self._run)
        self.close_btn = QPushButton("Close")
        self.close_btn.clicked.connect(self.close)
        buttons.addWidget(self.run_btn)
        buttons.addWidget(self.close_btn)

        layout = QVBoxLayout(self)
        layout.addLayout(grid)
        layout.addWidget(self.log, 1)
        layout.addLayout(buttons)

    # -- file pickers --------------------------------------------------------

    def _browse_input(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "ASCII forcefield grid", self.input_edit.text(),
            "ASCII grids (*.dat *.txt);;All files (*)")
        if path:
            self.input_edit.setText(path)
            if not self.output_edit.text():
                self.output_edit.setText(
                    str(Path(path).with_suffix("")) + ".v2.bin")

    def _browse_output(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Output binary", self.output_edit.text(),
            "FSMP forcefields (*.bin)")
        if path:
            self.output_edit.setText(path)

    def _browse_exe(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "pack executable", self.exe_edit.text(),
            "Executables (*.exe *.out);;All files (*)")
        if path:
            self.exe_edit.setText(path)

    # -- running -------------------------------------------------------------

    def _run(self) -> None:
        src = self.input_edit.text().strip()
        dst = self.output_edit.text().strip()
        exe = self.exe_edit.text().strip()
        if not Path(src).is_file():
            QMessageBox.warning(self, "No input", f"Input file not found:\n{src}")
            return
        if not dst:
            QMessageBox.warning(self, "No output", "Choose the output file.")
            return
        if not Path(exe).is_file():
            QMessageBox.warning(self, "No converter",
                                "Point to the pack executable (pack.exe from "
                                "the release zip, or build it with 'make pack').")
            return
        if Path(dst).exists():
            answer = QMessageBox.question(self, "Overwrite",
                                          f"{dst}\nalready exists. Overwrite?")
            if answer != QMessageBox.Yes:
                return
        QSettings().setValue("pack/exe", exe)

        args = []
        if self.float_box.isChecked():
            args.append("--float")
        if self.stride_spin.value() > 1:
            args += ["--stride", str(self.stride_spin.value())]
        args += [src, dst]
        if self.fold_spin.value() > 0:
            args.append(str(self.fold_spin.value()))

        self.success = False
        self.output_path = dst
        self.log.clear()
        self.log.appendPlainText(f"> {exe} {' '.join(args)}\n")
        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._append_output)
        self.process.finished.connect(self._finished)
        self.process.errorOccurred.connect(
            lambda err: self.log.appendPlainText(f"[process error: {err}]"))
        self.run_btn.setEnabled(False)
        self.close_btn.setText("Cancel")
        self.process.start(exe, args)

    def _append_output(self) -> None:
        text = bytes(self.process.readAllStandardOutput()).decode(errors="replace")
        self.log.appendPlainText(text.rstrip("\n"))

    def _finished(self, code, _status) -> None:
        self.log.appendPlainText(f"\n[finished with exit code {code}]")
        self.success = code == 0
        self.run_btn.setEnabled(True)
        self.close_btn.setText("Close")
        self.process = None

    def closeEvent(self, event) -> None:
        if self.process is not None and self.process.state() != QProcess.NotRunning:
            answer = QMessageBox.question(
                self, "Conversion running",
                "The conversion is still running. Stop it and close?")
            if answer != QMessageBox.Yes:
                event.ignore()
                return
            self.process.kill()
            self.process.waitForFinished(3000)
        event.accept()
