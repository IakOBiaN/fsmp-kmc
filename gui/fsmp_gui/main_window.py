"""Main window: start page <-> project view with the workflow tabs."""

from pathlib import Path

from PySide6.QtCore import QSettings, QStandardPaths, Qt
from PySide6.QtGui import QAction, QIcon
from PySide6.QtWidgets import (QDialog, QDialogButtonBox, QFileDialog,
                               QGridLayout, QHBoxLayout, QLabel, QLineEdit,
                               QMainWindow, QMessageBox, QPushButton,
                               QStackedWidget, QTabWidget, QVBoxLayout,
                               QWidget)

from . import theme
from .glyph import model_glyph
from .project import Project, ProjectError, safe_filename
from .start_page import ASSETS, StartPage
from .tabs.create_potential_tab import CreatePotentialTab
from .tabs.molecule_model_tab import MoleculeModelTab
from .tabs.potentials_tab import PotentialsTab
from .tabs.run_tab import RunTab
from .tabs.simulation_cell_tab import SimulationCellTab
from .tabs.unit_cell_tab import UnitCellTab

MAX_RECENT = 8


class NewProjectDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("New project")
        self.setMinimumWidth(520)

        grid = QGridLayout()
        grid.addWidget(QLabel("Name"), 0, 0)
        self.name_edit = QLineEdit()
        self.name_edit.setPlaceholderText("my-monolayer")
        grid.addWidget(self.name_edit, 0, 1, 1, 2)

        grid.addWidget(QLabel("Location"), 1, 0)
        documents = QStandardPaths.writableLocation(QStandardPaths.DocumentsLocation)
        self.location_edit = QLineEdit(documents)
        grid.addWidget(self.location_edit, 1, 1)
        browse = QPushButton("Browse…")
        browse.clicked.connect(self._browse)
        grid.addWidget(browse, 1, 2)

        self.preview = QLabel(" ")
        self.preview.setProperty("dim", True)

        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(grid)
        layout.addWidget(self.preview)
        layout.addWidget(self.buttons)

        self.name_edit.textChanged.connect(self._validate)
        self.location_edit.textChanged.connect(self._validate)
        self._validate()

    def _browse(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Project location",
                                                self.location_edit.text())
        if path:
            self.location_edit.setText(path)

    def target_path(self) -> Path | None:
        name = self.name_edit.text().strip()
        location = self.location_edit.text().strip()
        if not name or not location:
            return None
        return Path(location) / safe_filename(name)

    def _validate(self) -> None:
        target = self.target_path()
        ok = False
        if target is None:
            self.preview.setText("Enter a project name and location")
        elif target.exists() and any(target.iterdir()):
            self.preview.setText(f"Folder already exists and is not empty: {target}")
        elif not target.parent.is_dir():
            self.preview.setText(f"Location does not exist: {target.parent}")
        else:
            self.preview.setText(f"The project will be created in  {target}")
            ok = True
        self.buttons.button(QDialogButtonBox.Ok).setEnabled(ok)


class ProjectView(QWidget):
    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self.project = project

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 8, 12, 8)

        header = QHBoxLayout()
        name = QLabel(project.name)
        name.setProperty("gold", True)
        name.setStyleSheet(f"font-size: 13pt; font-weight: 600; color: {theme.GOLD};"
                           "background: transparent;")
        header.addWidget(name)
        path = QLabel(str(project.root))
        path.setProperty("dim", True)
        header.addWidget(path)
        header.addStretch(1)
        layout.addLayout(header)

        self.tabs = QTabWidget()
        self.model_tab = MoleculeModelTab(project)
        self.tabs.addTab(self.model_tab, "1  Molecule model")
        self.create_tab = CreatePotentialTab(project)
        self.tabs.addTab(self.create_tab, "2  Create potential")
        self.potentials_tab = PotentialsTab(project)
        self.tabs.addTab(self.potentials_tab, "3  Potentials")
        self.unit_cell_tab = UnitCellTab(project)
        self.tabs.addTab(self.unit_cell_tab, "4  Unit cell")
        self.sim_cell_tab = SimulationCellTab(project)
        self.tabs.addTab(self.sim_cell_tab, "5  Simulation cell")
        self.run_tab = RunTab(project)
        self.tabs.addTab(self.run_tab, "6  Run")
        layout.addWidget(self.tabs, 1)
        self.tabs.currentChanged.connect(self._tab_changed)
        # the project model can change on tab 1; downstream tabs must re-read it
        self.model_tab.projectModelChanged.connect(self.create_tab.refresh)
        self.model_tab.projectModelChanged.connect(
            lambda: self.unit_cell_tab.canvas.set_glyph(model_glyph(project)))
        # a freshly generated potential should show up on the Potentials tab
        self.create_tab.site_page.potentialGenerated.connect(
            self.potentials_tab.refresh)

    def _tab_changed(self, index: int) -> None:
        # tabs read shared project state, so refresh on activation
        widget = self.tabs.widget(index)
        if widget is self.create_tab:
            self.create_tab.refresh()
        elif widget is self.potentials_tab:
            self.potentials_tab.refresh()
        elif widget is self.unit_cell_tab:
            self.unit_cell_tab.refresh()
        elif widget is self.sim_cell_tab:
            self.sim_cell_tab.refresh()
        elif widget is self.run_tab:
            self.run_tab.refresh()


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("FSMP-kMC Studio")
        self.setWindowIcon(QIcon(str(ASSETS / "logo-mark.svg")))
        self.project: Project | None = None
        self.project_view: ProjectView | None = None

        self.stack = QStackedWidget()
        self.start_page = StartPage()
        self.start_page.newProjectRequested.connect(self.new_project)
        self.start_page.openProjectRequested.connect(self.open_project)
        self.start_page.recentProjectRequested.connect(self.open_project_at)
        self.stack.addWidget(self.start_page)
        self.setCentralWidget(self.stack)

        self._build_menu()
        self.statusBar().showMessage("Ready")
        self._refresh_recent()

    def _build_menu(self) -> None:
        file_menu = self.menuBar().addMenu("&File")
        for text, shortcut, slot in (
                ("&New project…", "Ctrl+N", self.new_project),
                ("&Open project…", "Ctrl+O", self.open_project)):
            action = QAction(text, self)
            action.setShortcut(shortcut)
            action.triggered.connect(slot)
            file_menu.addAction(action)
        file_menu.addSeparator()
        self.close_action = QAction("&Close project", self)
        self.close_action.setShortcut("Ctrl+W")
        self.close_action.setEnabled(False)
        self.close_action.triggered.connect(self.close_project)
        file_menu.addAction(self.close_action)
        file_menu.addSeparator()
        quit_action = QAction("E&xit", self)
        quit_action.setShortcut("Ctrl+Q")
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)

    # -- recent projects ---------------------------------------------------

    def _recent(self) -> list[str]:
        value = QSettings().value("recent_projects", [])
        if isinstance(value, str):
            value = [value]
        return list(value or [])

    def _add_recent(self, path: str) -> None:
        recent = [p for p in self._recent() if p != path]
        recent.insert(0, path)
        QSettings().setValue("recent_projects", recent[:MAX_RECENT])
        self._refresh_recent()

    def _refresh_recent(self) -> None:
        self.start_page.set_recent(self._recent())

    # -- project lifecycle -------------------------------------------------

    def new_project(self) -> None:
        dialog = NewProjectDialog(self)
        if dialog.exec() != QDialog.Accepted:
            return
        target = dialog.target_path()
        try:
            project = Project.create(target, dialog.name_edit.text().strip())
        except (ProjectError, OSError) as e:
            QMessageBox.warning(self, "Cannot create project", str(e))
            return
        self._show_project(project)

    def open_project(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Open project folder")
        if path:
            self.open_project_at(path)

    def open_project_at(self, path: str) -> None:
        try:
            project = Project.open(path)
        except (ProjectError, OSError) as e:
            QMessageBox.warning(self, "Cannot open project", str(e))
            return
        self._show_project(project)

    def _show_project(self, project: Project) -> None:
        if self.project_view is not None:
            self.stack.removeWidget(self.project_view)
            self.project_view.deleteLater()
        self.project = project
        self.project_view = ProjectView(project)
        self.project_view.model_tab.statusMessage.connect(
            self.statusBar().showMessage)
        self.project_view.create_tab.statusMessage.connect(
            self.statusBar().showMessage)
        self.project_view.potentials_tab.statusMessage.connect(
            self.statusBar().showMessage)
        self.project_view.unit_cell_tab.statusMessage.connect(
            self.statusBar().showMessage)
        self.project_view.sim_cell_tab.statusMessage.connect(
            self.statusBar().showMessage)
        self.project_view.run_tab.statusMessage.connect(
            self.statusBar().showMessage)
        self.stack.addWidget(self.project_view)
        self.stack.setCurrentWidget(self.project_view)
        self.close_action.setEnabled(True)
        self.setWindowTitle(f"{project.name} — FSMP-kMC Studio")
        self._add_recent(str(project.root))
        self.statusBar().showMessage(f"Project: {project.root}")

    def close_project(self) -> None:
        if self.project_view is not None:
            self.stack.removeWidget(self.project_view)
            self.project_view.deleteLater()
            self.project_view = None
        self.project = None
        self.close_action.setEnabled(False)
        self.setWindowTitle("FSMP-kMC Studio")
        self.stack.setCurrentWidget(self.start_page)
        self._refresh_recent()
