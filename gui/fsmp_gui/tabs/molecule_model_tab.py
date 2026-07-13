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
        # the atomistic model is the entry point: every visualization needs
        # it, so it is the first thing a fresh project asks for
        self.subtabs.addTab(self.atomistic_tab, "Atomistic model")
        self.subtabs.addTab(self.site_tab, "Site model")
        layout.addWidget(self.subtabs)

        for tab in (self.atomistic_tab, self.site_tab):
            tab.statusMessage.connect(self.statusMessage)
            tab.projectModelChanged.connect(self._model_changed)
        self.refresh_gating()

    def refresh_gating(self) -> None:
        """The site subtab opens up once the atomistic model is attached."""
        self.subtabs.setTabEnabled(self.subtabs.indexOf(self.site_tab),
                                   self.project.atomistic is not None)

    def _model_changed(self) -> None:
        # the two subtabs share one project slot; keep both panels in sync
        self.atomistic_tab.refresh_project_model()
        self.site_tab.refresh_project_model()
        self.refresh_gating()
        self.projectModelChanged.emit()
