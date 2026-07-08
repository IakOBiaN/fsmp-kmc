"""Tab 1 container: the project's molecule model, either atomistic geometry
or a coarse-grained site model, on two subtabs sharing one project slot."""

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QTabWidget, QVBoxLayout, QWidget

from ..project import Project
from .molecule_tab import MoleculeTab
from .site_tab import SiteTab


class MoleculeModelTab(QWidget):
    statusMessage = Signal(str)
    projectModelChanged = Signal()

    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self.project = project

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.subtabs = QTabWidget()
        self.atomistic_tab = MoleculeTab(project)
        self.site_tab = SiteTab(project)
        self.subtabs.addTab(self.atomistic_tab, "Atomistic model")
        self.subtabs.addTab(self.site_tab, "Site model")
        self.subtabs.setCurrentWidget(self.site_tab)  # charge models are the focus
        layout.addWidget(self.subtabs)

        for tab in (self.atomistic_tab, self.site_tab):
            tab.statusMessage.connect(self.statusMessage)
            tab.projectModelChanged.connect(self._model_changed)

    def _model_changed(self) -> None:
        # the two subtabs share one project slot; keep both panels in sync
        self.atomistic_tab.refresh_project_model()
        self.site_tab.refresh_project_model()
        self.projectModelChanged.emit()
