"""Tab 6 "Run": launch and follow production runs.

Left: the run form (the T and u_m loops and the Monte Carlo settings; the
potential, the unit cell and the elongated-cell parameters come from the
project). Start creates <project>/runs/<name>/ and launches the engine
detached, so runs survive closing the GUI.

Right: one card per run with live progress (current point, T, u_m, percent),
and for the selected run the details: statistics plots with selectable
columns, the lazy trajectory viewer and the engine log.
"""

from datetime import datetime
from pathlib import Path

from PySide6.QtCore import Qt, QTimer, QUrl, Signal
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFrame,
                               QGridLayout, QGroupBox, QHBoxLayout, QLabel,
                               QLineEdit, QMessageBox, QPlainTextEdit,
                               QProgressBar, QPushButton, QScrollArea,
                               QSpinBox, QSplitter, QTabWidget, QVBoxLayout,
                               QWidget)

from .. import runs, theme
from ..engine import find_engine
from ..plot import PlotWidget
from ..project import Project
from ..sitemodel import SiteModel
from ..trajectory import TrajectoryViewer

POLL_MS = 1000
LOG_CAP = 2 * 1024 * 1024   # show at most the last 2 MB of a log

STATE_COLORS = {runs.RUNNING: theme.GOLD, runs.DONE: theme.ACCENT,
                runs.FAILED: theme.DANGER, runs.STOPPED: theme.TEXT_DIM,
                runs.INTERRUPTED: theme.DANGER}


class RunCard(QFrame):
    """One run in the list: label, state, live progress."""

    selected = Signal(object)   # self

    def __init__(self, run_dir: Path, parent=None):
        super().__init__(parent)
        self.run_dir = run_dir
        self.watch = runs.LogWatch(run_dir / runs.LOG)
        self.alive: bool | None = None   # verified on scan only
        try:
            params = runs.read_params(run_dir / runs.PARAMS)
            self.ts, self.ums = runs.loop_points(params)
        except (OSError, ValueError, KeyError):
            self.ts, self.ums = [0.0], [0.0]
        self.total = len(self.ts) * len(self.ums)

        self.setProperty("panel", True)
        self.setFrameShape(QFrame.StyledPanel)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 8, 10, 8)
        layout.setSpacing(4)

        top = QHBoxLayout()
        name = QLabel(f"<b>{runs.run_label(run_dir)}</b>")
        top.addWidget(name)
        self.state_label = QLabel()
        top.addWidget(self.state_label)
        top.addStretch(1)
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.clicked.connect(self._stop)
        top.addWidget(self.stop_btn)
        folder = QPushButton("Folder")
        folder.clicked.connect(lambda: QDesktopServices.openUrl(
            QUrl.fromLocalFile(str(self.run_dir))))
        top.addWidget(folder)
        layout.addLayout(top)

        self.bar = QProgressBar()
        self.bar.setRange(0, 1000)
        self.bar.setTextVisible(False)
        self.bar.setFixedHeight(8)
        layout.addWidget(self.bar)
        self.status = QLabel(" ")
        self.status.setProperty("dim", True)
        layout.addWidget(self.status)
        self.poll()

    def mousePressEvent(self, event):
        self.selected.emit(self)
        super().mousePressEvent(event)

    def set_highlight(self, on: bool) -> None:
        color = theme.ACCENT if on else theme.BORDER
        self.setStyleSheet(f"QFrame[panel='true'] {{ border: 1px solid "
                           f"{color}; border-radius: 6px; "
                           f"background: {theme.BG_RAISED}; }}")

    # -- state ---------------------------------------------------------------

    def state(self) -> str:
        return self.watch.state(self.alive)

    def poll(self) -> None:
        changed = self.watch.poll()
        if changed or self.state_label.text() == "":
            self._render()

    def verify_alive(self) -> None:
        """Process check, done on scan only."""
        if self.watch.exit_code is None:
            self.alive = runs.is_alive(self.run_dir)
            self._render()

    def _render(self) -> None:
        state = self.state()
        self.state_label.setText(
            f"<span style='color:{STATE_COLORS[state]}'>{state}</span>")
        self.stop_btn.setVisible(state == runs.RUNNING)
        w = self.watch
        if state == runs.RUNNING and w.points_started == 0:
            self.bar.setValue(0)
            self.status.setText("optimizing the unit cell…")
            return
        if state == runs.RUNNING:
            index = w.points_started - 1
            t, um = runs.point_at(self.ts, self.ums, min(index, self.total - 1))
            self.bar.setValue(int((index + w.percent / 100.0)
                                  / self.total * 1000))
            self.status.setText(
                f"point {w.points_started}/{self.total}   ·   T = {t:g} K   ·"
                f"   u_m = {um / 1000.0:g} kJ/mol   ·   {w.percent} %")
            return
        self.bar.setValue(1000 if state == runs.DONE else
                          int(max(w.points_started - 1, 0) / self.total * 1000))
        done = f"{max(w.points_started + (-1 if state != runs.DONE else 0), 0)}"
        if state == runs.DONE:
            self.status.setText(f"{self.total} points finished")
        elif state == runs.FAILED:
            self.status.setText(f"exit code {w.exit_code}; see the log")
        else:
            self.status.setText(f"stopped after {done} of {self.total} points")

    def _stop(self) -> None:
        if QMessageBox.question(self, "Stop run",
                                f"Stop '{runs.run_label(self.run_dir)}'?"
                                ) != QMessageBox.Yes:
            return
        runs.stop(self.run_dir)


class RunDetail(QWidget):
    """Plots, trajectory and log of the selected run."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.run_dir: Path | None = None
        self._stats_mtime = 0.0
        self._log_pos = 0

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)

        plots = QWidget()
        pv = QVBoxLayout(plots)
        row = QHBoxLayout()
        row.addWidget(QLabel("X"))
        self.x_box = QComboBox()
        row.addWidget(self.x_box, 1)
        row.addWidget(QLabel("Y"))
        self.y_box = QComboBox()
        row.addWidget(self.y_box, 1)
        pv.addLayout(row)
        self.plot = PlotWidget()
        pv.addWidget(self.plot, 1)
        self.x_box.currentIndexChanged.connect(self._replot)
        self.y_box.currentIndexChanged.connect(self._replot)
        self.tabs.addTab(plots, "Plots")

        self.viewer = TrajectoryViewer()
        self.tabs.addTab(self.viewer, "Trajectory")

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(20000)
        self.tabs.addTab(self.log, "Log")

        self._names: list = []
        self._rows: list = []

    def set_run(self, run_dir: Path | None) -> None:
        if run_dir == self.run_dir:
            return
        self.run_dir = run_dir
        self._stats_mtime = 0.0
        self._log_pos = 0
        self._names, self._rows = [], []
        self.log.clear()
        self.x_box.blockSignals(True)
        self.y_box.blockSignals(True)
        self.x_box.clear()
        self.y_box.clear()
        self.x_box.blockSignals(False)
        self.y_box.blockSignals(False)
        self.plot.set_data([], [], "", "")
        if run_dir is not None:
            self.viewer.set_path(run_dir / runs.TRAJECTORY)
            self.poll()

    def poll(self) -> None:
        """Pick up new statistics rows and log text."""
        if self.run_dir is None:
            return
        stats = self.run_dir / runs.STATISTICS
        try:
            mtime = stats.stat().st_mtime
        except OSError:
            mtime = 0.0
        if mtime != self._stats_mtime:
            self._stats_mtime = mtime
            self._names, self._rows = runs.read_statistics(stats)
            self._fill_boxes()
            self._replot()
        self._poll_log()
        if self.viewer.isVisible():
            self.viewer.poll()   # keep showing the newest frame

    def _fill_boxes(self) -> None:
        if not self._names or self.x_box.count() == len(self._names):
            return
        def pick(box, needle, fallback):
            for i, name in enumerate(self._names):
                if needle in name:
                    box.setCurrentIndex(i)
                    return
            box.setCurrentIndex(fallback)
        for box in (self.x_box, self.y_box):
            box.blockSignals(True)
            box.clear()
            box.addItems(self._names)
        # the house convention: chemical potential on x, pressure on y
        pick(self.x_box, "kMC, gas phase", 0)
        pick(self.y_box, "Analytical pressure", min(7, len(self._names) - 1))
        for box in (self.x_box, self.y_box):
            box.blockSignals(False)

    def _replot(self) -> None:
        xi, yi = self.x_box.currentIndex(), self.y_box.currentIndex()
        if not self._rows or xi < 0 or yi < 0:
            self.plot.set_data([], [], "", "")
            return
        self.plot.set_data([r[xi] for r in self._rows],
                           [r[yi] for r in self._rows],
                           self._names[xi], self._names[yi])

    def _poll_log(self) -> None:
        path = self.run_dir / runs.LOG
        try:
            size = path.stat().st_size
            if size <= self._log_pos:
                return
            with open(path, "rb") as f:
                if self._log_pos == 0 and size > LOG_CAP:
                    f.seek(size - LOG_CAP)
                else:
                    f.seek(self._log_pos)
                text = f.read().decode("utf-8", "replace")
                self._log_pos = f.tell()
        except OSError:
            return
        cursor = self.log.textCursor()
        cursor.movePosition(cursor.MoveOperation.End)
        cursor.insertText(text)
        self.log.setTextCursor(cursor)


class RunTab(QWidget):
    statusMessage = Signal(str)

    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self.project = project
        self.cards: dict[Path, RunCard] = {}
        self.current: RunCard | None = None

        root = QHBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        splitter = QSplitter(Qt.Horizontal)
        root.addWidget(splitter)
        splitter.addWidget(self._build_form())
        splitter.addWidget(self._build_right())
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([420, 880])

        self.timer = QTimer(self)
        self.timer.setInterval(POLL_MS)
        self.timer.timeout.connect(self._poll)

    # -- form ------------------------------------------------------------------

    def _dspin(self, lo, hi, val, step, dec=1, suffix=""):
        s = QDoubleSpinBox()
        s.setRange(lo, hi)
        s.setDecimals(dec)
        s.setSingleStep(step)
        s.setValue(val)
        if suffix:
            s.setSuffix(suffix)
        return s

    def _ispin(self, lo, hi, val):
        s = QSpinBox()
        s.setRange(lo, hi)
        s.setValue(val)
        return s

    def _build_form(self) -> QWidget:
        wrap = QScrollArea()
        wrap.setWidgetResizable(True)
        wrap.setFrameShape(QScrollArea.NoFrame)
        panel = QWidget()
        wrap.setWidget(panel)
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(0, 0, 8, 0)

        self.prereq = QLabel()
        self.prereq.setWordWrap(True)
        self.prereq.setStyleSheet(f"color: {theme.GOLD};")
        layout.addWidget(self.prereq)

        self.engine_label = QLabel()
        self.engine_label.setWordWrap(True)
        self.engine_label.setProperty("dim", True)
        layout.addWidget(self.engine_label)

        name_row = QHBoxLayout()
        name_row.addWidget(QLabel("Run name"))
        self.name_edit = QLineEdit()
        name_row.addWidget(self.name_edit, 1)
        layout.addLayout(name_row)

        loops = QGroupBox("T and u_m loops  (u_m is the outer loop)")
        grid = QGridLayout(loops)
        self.temp_from = self._dspin(1, 3000, 300, 10, suffix=" K")
        self.temp_to = self._dspin(1, 3000, 300, 10, suffix=" K")
        self.temp_step = self._dspin(0.1, 1000, 10, 1)
        self.um_from = self._dspin(-1e6, 1e6, 0, 1000, dec=0, suffix=" J/mol")
        self.um_to = self._dspin(-1e6, 1e6, -20000, 1000, dec=0,
                                 suffix=" J/mol")
        self.um_step = self._dspin(1, 1e6, 1000, 500, dec=0, suffix=" J/mol")
        for r, (text, w) in enumerate((("T from", self.temp_from),
                                       ("T to", self.temp_to),
                                       ("T step", self.temp_step),
                                       ("u_m from", self.um_from),
                                       ("u_m to", self.um_to),
                                       ("u_m step", self.um_step))):
            grid.addWidget(QLabel(text), r, 0)
            grid.addWidget(w, r, 1)
        self.points_label = QLabel()
        self.points_label.setProperty("dim", True)
        grid.addWidget(self.points_label, 6, 0, 1, 2)
        for w in (self.temp_from, self.temp_to, self.temp_step,
                  self.um_from, self.um_to, self.um_step):
            w.valueChanged.connect(self._update_points)
        layout.addWidget(loops)

        mc = QGroupBox("Monte Carlo")
        grid2 = QGridLayout(mc)
        self.nsteps = self._ispin(1, 1000000000, 500000)
        self.nsteps_eq = self._ispin(1, 1000000000, 250000)
        self.delta = self._dspin(0.01, 50, 2.0, 0.1, dec=2, suffix=" Å")
        self.delta_angle = self._dspin(0.1, 180, 60.0, 5, suffix=" °")
        self.kmc = QCheckBox("kinetic MC (kMC)")
        self.kmc.setChecked(True)
        self.widom = QCheckBox("Widom test")
        self.seed = self._ispin(0, 2147483647, 0)
        self.seed.setSpecialValueText("random")
        rows = (("MC steps (nSteps)", self.nsteps),
                ("relaxation (nStepsEq)", self.nsteps_eq),
                ("max shift", self.delta),
                ("max rotation", self.delta_angle),
                ("seed", self.seed))
        for r, (text, w) in enumerate(rows):
            grid2.addWidget(QLabel(text), r, 0)
            grid2.addWidget(w, r, 1)
        grid2.addWidget(self.kmc, len(rows), 0)
        grid2.addWidget(self.widom, len(rows), 1)
        layout.addWidget(mc)

        extra = QGroupBox("Reference area and pressure")
        grid3 = QGridLayout(extra)
        self.sigma_mode = QComboBox()
        self.sigma_mode.addItems(["manual", "min_dist", "molecule_area"])
        self.sigma = self._dspin(0.1, 100, 11.052, 0.1, dec=3, suffix=" Å")
        self.sigma_mode.currentTextChanged.connect(
            lambda mode: self.sigma.setEnabled(mode == "manual"))
        grid3.addWidget(QLabel("sigma_mode"), 0, 0)
        grid3.addWidget(self.sigma_mode, 0, 1)
        grid3.addWidget(QLabel("sigma"), 1, 0)
        grid3.addWidget(self.sigma, 1, 1)
        self.const_p = QCheckBox("constant pressure")
        self.const_p_value = self._dspin(-1e6, 1e6, 0, 1, dec=2, suffix=" mN/m")
        self.const_p_value.setEnabled(False)
        self.const_p.toggled.connect(self.const_p_value.setEnabled)
        grid3.addWidget(self.const_p, 2, 0)
        grid3.addWidget(self.const_p_value, 2, 1)
        layout.addWidget(extra)

        mask_box = QGroupBox("Stabilization mask (porous phases)")
        grid4 = QGridLayout(mask_box)
        self.mask_enable = QCheckBox("use the mask in this run")
        self.mask_enable.setToolTip("A lattice of free wells built from the "
                                    "initial structure that keeps a "
                                    "metastable polymorph intact")
        grid4.addWidget(self.mask_enable, 0, 0, 1, 2)
        self.mask_radius = self._dspin(0.1, 50, 3.0, 0.5, suffix=" Å")
        self.mask_ramp = self._dspin(0.1, 50, 2.0, 0.5, suffix=" Å")
        self.mask_penalty = self._dspin(1, 1e6, 25000, 1000, dec=0,
                                        suffix=" J/mol")
        for r, (text, w) in enumerate((("free radius", self.mask_radius),
                                       ("ramp width", self.mask_ramp),
                                       ("penalty", self.mask_penalty))):
            grid4.addWidget(QLabel(text), r + 1, 0)
            grid4.addWidget(w, r + 1, 1)
            w.setEnabled(False)
            self.mask_enable.toggled.connect(w.setEnabled)
        layout.addWidget(mask_box)

        start = QPushButton("Start run")
        start.setProperty("primary", True)
        start.clicked.connect(self._start)
        layout.addWidget(start)
        layout.addStretch(1)
        self._update_points()
        return wrap

    def _update_points(self) -> None:
        ts = runs.loop_values(self.temp_from.value(), self.temp_to.value(),
                              self.temp_step.value())
        ums = runs.loop_values(self.um_from.value(), self.um_to.value(),
                               self.um_step.value())
        self.points_label.setText(
            f"{len(ts)} temperatures × {len(ums)} u_m values = "
            f"{len(ts) * len(ums)} points")

    # -- the run list and details ------------------------------------------------

    def _build_right(self) -> QWidget:
        splitter = QSplitter(Qt.Vertical)

        wrap = QScrollArea()
        wrap.setWidgetResizable(True)
        wrap.setFrameShape(QScrollArea.NoFrame)
        inner = QWidget()
        self.cards_layout = QVBoxLayout(inner)
        self.cards_layout.setContentsMargins(0, 0, 0, 0)
        self.cards_layout.setSpacing(6)
        self.empty_label = QLabel("No runs yet. Fill the form and press "
                                  "'Start run'.")
        self.empty_label.setProperty("dim", True)
        self.empty_label.setAlignment(Qt.AlignCenter)
        self.cards_layout.addWidget(self.empty_label)
        self.cards_layout.addStretch(1)
        wrap.setWidget(inner)
        splitter.addWidget(wrap)

        self.detail = RunDetail()
        splitter.addWidget(self.detail)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([260, 480])
        return splitter

    # -- lifecycle -----------------------------------------------------------------

    def refresh(self) -> None:
        """Tab activation: scan run folders, load form defaults, start polling."""
        self._load_form_defaults()
        self._check_prereqs()
        for run_dir in runs.list_runs(self.project):
            if run_dir not in self.cards:
                self._add_card(run_dir, verify=True)
        if self.name_edit.text() == "":
            self.name_edit.setPlaceholderText(
                datetime.now().strftime("run-%Y%m%d-%H%M%S"))
        self.timer.start()

    def hideEvent(self, event):
        self.timer.stop()
        super().hideEvent(event)

    def showEvent(self, event):
        super().showEvent(event)
        if self.cards:
            self.timer.start()

    def _check_prereqs(self) -> None:
        missing = []
        if self.project.potential is None:
            missing.append("potential (tab 3)")
        if self.project.unit_cell is None:
            missing.append("unit cell (tab 4)")
        if self.project.simulation_cell is None:
            missing.append("simulation cell (tab 5)")
        if self.project.atomistic is None:
            missing.append("atomistic model (tab 1)")
        command = find_engine()
        if command is None:
            missing.append("engine (a release fsmp.exe in the repository "
                           "root, or make)")
            self.engine_label.setText("")
        else:
            self.engine_label.setText(f"Engine: {command[0]}")
        self.engine_label.setVisible(bool(self.engine_label.text()))
        self.prereq.setText("Missing before a run can start: "
                            + ", ".join(missing) if missing else "")
        self.prereq.setVisible(bool(missing))

    def _load_form_defaults(self) -> None:
        form = self.project.simulation
        if form is None:
            # sigma default: the LJ diameter of the site model when present
            entry = self.project.site
            if entry is not None:
                try:
                    sm = SiteModel.load(self.project.model_path(entry))
                    sigma = max((s.sigma for s in sm.sites if s.is_lj),
                                default=0.0)
                    if sigma > 0:
                        self.sigma.setValue(sigma)
                except (OSError, ValueError):
                    pass
            return
        widgets = {
            "temp_from": self.temp_from, "temp_to": self.temp_to,
            "temp_step": self.temp_step, "um_from": self.um_from,
            "um_to": self.um_to, "um_step": self.um_step,
            "nSteps": self.nsteps, "nStepsEq": self.nsteps_eq,
            "delta": self.delta, "delta_angle": self.delta_angle,
            "sigma": self.sigma, "seed": self.seed,
            "constant_pressure_value": self.const_p_value,
            "mask_free_radius": self.mask_radius,
            "mask_ramp_width": self.mask_ramp,
            "mask_penalty": self.mask_penalty,
        }
        for key, w in widgets.items():
            if key in form:
                w.setValue(form[key])
        self.kmc.setChecked(bool(form.get("kMC", True)))
        self.widom.setChecked(bool(form.get("widom", False)))
        self.const_p.setChecked(bool(form.get("constant_pressure", False)))
        self.mask_enable.setChecked(bool(form.get("mask", False)))
        mode = form.get("sigma_mode", "manual")
        if mode in ("manual", "min_dist", "molecule_area"):
            self.sigma_mode.setCurrentText(mode)

    def _form(self) -> dict:
        return {
            "temp_from": self.temp_from.value(),
            "temp_to": self.temp_to.value(),
            "temp_step": self.temp_step.value(),
            "um_from": self.um_from.value(),
            "um_to": self.um_to.value(),
            "um_step": self.um_step.value(),
            "nSteps": self.nsteps.value(),
            "nStepsEq": self.nsteps_eq.value(),
            "kMC": self.kmc.isChecked(),
            "delta": self.delta.value(),
            "delta_angle": self.delta_angle.value(),
            "widom": self.widom.isChecked(),
            "sigma_mode": self.sigma_mode.currentText(),
            "sigma": self.sigma.value(),
            "constant_pressure": self.const_p.isChecked(),
            "constant_pressure_value": self.const_p_value.value(),
            "mask": self.mask_enable.isChecked(),
            "mask_free_radius": self.mask_radius.value(),
            "mask_ramp_width": self.mask_ramp.value(),
            "mask_penalty": self.mask_penalty.value(),
            "seed": self.seed.value(),
        }

    # -- actions ----------------------------------------------------------------------

    def _start(self) -> None:
        if self.nsteps_eq.value() >= self.nsteps.value():
            QMessageBox.warning(self, "Bad steps",
                                "nStepsEq must be below nSteps.")
            return
        label = self.name_edit.text().strip() or self.name_edit.placeholderText()
        form = self._form()
        try:
            run_dir = runs.create_run(self.project, label, form)
        except (runs.RunError, OSError) as e:
            QMessageBox.warning(self, "Cannot start the run", str(e))
            return
        self.project.set_simulation(form)
        self.name_edit.clear()
        self.name_edit.setPlaceholderText(
            datetime.now().strftime("run-%Y%m%d-%H%M%S"))
        card = self._add_card(run_dir, verify=False, to_top=True)
        card.alive = True
        self._select(card)
        self.timer.start()
        self.statusMessage.emit(f"Run started in {run_dir}")

    def _add_card(self, run_dir: Path, verify: bool,
                  to_top: bool = False) -> RunCard:
        card = RunCard(run_dir)
        card.selected.connect(self._select)
        self.cards[run_dir] = card
        self.empty_label.setVisible(False)
        index = 0 if to_top else self.cards_layout.count() - 1
        self.cards_layout.insertWidget(index, card)
        card.set_highlight(False)
        if verify:
            card.verify_alive()
        if self.current is None:
            self._select(card)
        return card

    def _select(self, card: RunCard) -> None:
        if self.current is card:
            return
        if self.current is not None:
            self.current.set_highlight(False)
        self.current = card
        card.set_highlight(True)
        self.detail.set_run(card.run_dir)

    def _poll(self) -> None:
        for card in self.cards.values():
            card.poll()
        if self.current is not None:
            self.detail.poll()
