"""Site model editor (subtab of "Molecule model").

Place Lennard-Jones centers and point charges as circles, edit their
parameters, and attach the model to the project."""

from pathlib import Path

from PySide6.QtCore import QSettings, Qt, Signal
from PySide6.QtWidgets import (QAbstractItemView, QButtonGroup, QComboBox,
                               QFileDialog, QFrame, QHBoxLayout, QHeaderView,
                               QInputDialog, QLabel, QMessageBox, QPushButton,
                               QSplitter, QTableView, QToolButton, QVBoxLayout,
                               QWidget)

from ..canvas import Mode, SiteCanvas
from ..project import Project
from ..site_table import SiteTableModel
from ..sitemodel import SUFFIX, Site, SiteModel

# add-tool templates: label -> (default charge, is LJ)
TEMPLATES = {
    "LJ center": ("LJ", 0.0, True),
    "+ charge": ("q+", 1.0, False),
    "− charge": ("q-", -1.0, False),
}
_DEFAULT_LJ = (3.4, 50.0)  # sigma (Å), epsilon (K) for a fresh LJ center


class SiteTab(QWidget):
    statusMessage = Signal(str)
    projectModelChanged = Signal()

    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self.project = project
        self.model = SiteTableModel(self)
        self.model.changed.connect(self._on_model_changed)
        self._dirty = False
        self._current_name: str | None = None

        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.addLayout(self._build_toolbar())

        splitter = QSplitter(Qt.Horizontal)
        root.addWidget(splitter, 1)

        self.canvas = SiteCanvas(self.model)
        self.canvas.cursorMoved.connect(
            lambda x, y: self.statusMessage.emit(f"x = {x:.2f} Å,  y = {y:.2f} Å"))
        self.canvas.addRequested.connect(self._add_site_from_canvas)
        splitter.addWidget(self.canvas)
        splitter.addWidget(self._build_side_panel())
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 0)
        splitter.setSizes([880, 400])

        self._set_mode(Mode.SELECT)
        self.refresh_project_model()
        self._update_info()

    # -- construction ------------------------------------------------------

    def _build_toolbar(self) -> QHBoxLayout:
        bar = QHBoxLayout()
        bar.setSpacing(6)

        for text, slot in (("New", self.new_model),
                           ("Open…", self.open_model),
                           ("Save as…", self.save_model)):
            b = QPushButton(text)
            b.clicked.connect(slot)
            bar.addWidget(b)

        bar.addSpacing(16)

        self.mode_group = QButtonGroup(self)
        self.mode_buttons = {}
        for mode, text, tip in ((Mode.SELECT, "Select", "Select and drag sites"),
                                (Mode.ADD, "Add site", "Click empty space to add a site"),
                                (Mode.DELETE, "Delete", "Click a site to delete it")):
            b = QToolButton()
            b.setText(text)
            b.setToolTip(tip)
            b.setCheckable(True)
            b.clicked.connect(lambda _=False, m=mode: self._set_mode(m))
            self.mode_group.addButton(b)
            self.mode_buttons[mode] = b
            bar.addWidget(b)

        self.type_box = QComboBox()
        self.type_box.addItems(TEMPLATES.keys())
        self.type_box.setToolTip("What to place in 'Add site' mode")
        bar.addWidget(self.type_box)

        bar.addSpacing(16)
        reset = QPushButton("Fit view")
        reset.clicked.connect(lambda: self.canvas.reset_view())
        bar.addWidget(reset)

        bar.addStretch(1)

        attach = QPushButton("Use in project")
        attach.setProperty("primary", True)
        attach.setToolTip("Save this site model as the project model")
        attach.clicked.connect(self.use_in_project)
        bar.addWidget(attach)
        return bar

    def _build_side_panel(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(8, 0, 0, 0)

        caption = QLabel("Sites")
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
        caption2 = QLabel("Project model")
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
        self.slot_open.clicked.connect(self._open_project_model)
        self.slot_clear = QPushButton("Remove")
        self.slot_clear.clicked.connect(self._clear_project_model)
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
        sm = self.model.site_model
        text = sm.summary()
        if sm.sites:
            text += f"   ·   total charge {sm.total_charge():+.4f} e"
        self.info.setText(text)

    def _confirm_discard(self) -> bool:
        if not self._dirty or not self.model.site_model.sites:
            return True
        answer = QMessageBox.question(
            self, "Unsaved changes",
            "The current site model has unsaved changes. Discard them?")
        return answer == QMessageBox.Yes

    def _set_mode(self, mode: Mode) -> None:
        self.mode_buttons[mode].setChecked(True)
        self.canvas.set_mode(mode)

    # -- editing actions ---------------------------------------------------

    def _template_site(self, x: float, y: float) -> Site:
        label, q, is_lj = TEMPLATES[self.type_box.currentText()]
        if is_lj:
            sigma, eps = _DEFAULT_LJ
            return Site(label, x, y, 0.0, 0.0, eps, sigma)
        return Site(label, x, y, 0.0, q, 0.0, 0.0)

    def _add_site_from_canvas(self, x: float, y: float) -> None:
        self.model.add_site(self._template_site(x, y))

    def _add_row(self) -> None:
        self.model.add_site(self._template_site(0.0, 0.0))

    def _remove_selected_rows(self) -> None:
        rows = [i.row() for i in self.table.selectionModel().selectedRows()]
        self.model.remove_rows(rows)

    # -- file actions ------------------------------------------------------

    def _last_dir(self) -> str:
        return QSettings().value("site/last_dir", str(Path.home()))

    def _remember_dir(self, path: str) -> None:
        QSettings().setValue("site/last_dir", str(Path(path).parent))

    def new_model(self) -> None:
        if not self._confirm_discard():
            return
        self.model.set_site_model(SiteModel())
        self._current_name = None
        self._dirty = False
        self.statusMessage.emit("New site model")

    def open_model(self) -> None:
        if not self._confirm_discard():
            return
        path, _ = QFileDialog.getOpenFileName(
            self, "Open site model", self._last_dir(),
            f"Site models (*{SUFFIX});;All files (*)")
        if not path:
            return
        self._load_file(path)

    def _load_file(self, path: str, name: str | None = None) -> None:
        try:
            sm = SiteModel.load(path)
        except (OSError, ValueError) as e:
            QMessageBox.warning(self, "Cannot open site model", f"{path}\n\n{e}")
            return
        self.model.set_site_model(sm)
        self._current_name = name if name is not None else Path(path).stem
        self._dirty = False
        self._remember_dir(path)
        self.canvas.reset_view()
        self.statusMessage.emit(f"Opened {path}")

    def save_model(self) -> None:
        if not self.model.site_model.sites:
            QMessageBox.information(self, "Nothing to save", "The model is empty.")
            return
        suggestion = str(Path(self._last_dir())
                         / f"{self._current_name or 'sites'}{SUFFIX}")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save site model", suggestion, f"Site models (*{SUFFIX})")
        if not path:
            return
        try:
            self.model.site_model.save(path)
        except OSError as e:
            QMessageBox.warning(self, "Cannot save site model", str(e))
            return
        self._current_name = Path(path).stem
        self._dirty = False
        self._remember_dir(path)
        self.statusMessage.emit(f"Saved {path}")

    # -- project integration -----------------------------------------------

    def refresh_project_model(self) -> None:
        entry = self.project.model
        if entry is None:
            self.slot_label.setText("Not set yet. Build a site model and press "
                                    "'Use in project'.")
            self.slot_open.setEnabled(False)
            self.slot_clear.setEnabled(False)
        elif entry["kind"] == "site":
            self.slot_label.setText(f"{entry['name']}   ·   site   ·   "
                                    f"{entry['file']}")
            self.slot_open.setEnabled(True)
            self.slot_clear.setEnabled(True)
        else:
            self.slot_label.setText(f"{entry['name']}   ·   atomistic model. "
                                    "Using a site model here will replace it.")
            self.slot_open.setEnabled(False)
            self.slot_clear.setEnabled(True)

    def use_in_project(self) -> None:
        sm = self.model.site_model
        if not sm.sites:
            QMessageBox.information(self, "Empty model",
                                    "Add some sites before using the model in "
                                    "the project.")
            return
        default = self._current_name or "sites"
        name, ok = QInputDialog.getText(self, "Use site model in project",
                                        "Model name:", text=default)
        if not ok or not name.strip():
            return
        name = name.strip()
        old = self.project.model
        if old is not None:
            kind = "atomistic model" if old["kind"] == "atomistic" else "site model"
            answer = QMessageBox.question(
                self, "Replace project model",
                f"The project {kind} '{old['name']}' will be replaced. Continue?")
            if answer != QMessageBox.Yes:
                return
        self.project.set_site_model(name, sm)
        self._current_name = name
        self._dirty = False
        self.refresh_project_model()
        self.projectModelChanged.emit()
        self.statusMessage.emit(f"'{name}' is now the project model")

    def _open_project_model(self) -> None:
        entry = self.project.model
        if entry is None or entry["kind"] != "site":
            return
        if not self._confirm_discard():
            return
        self._load_file(str(self.project.model_path(entry)), name=entry["name"])

    def _clear_project_model(self) -> None:
        entry = self.project.model
        if entry is None:
            return
        answer = QMessageBox.question(
            self, "Remove project model",
            f"Remove '{entry['name']}' from the project?\n"
            "The file inside the project folder will be deleted.")
        if answer != QMessageBox.Yes:
            return
        self.project.clear_model()
        self.refresh_project_model()
        self.projectModelChanged.emit()
        self.statusMessage.emit("Project model removed")
