"""Start page: logo, new/open project, recent projects."""

from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtSvgWidgets import QSvgWidget
from PySide6.QtWidgets import (QHBoxLayout, QLabel, QListWidget,
                               QListWidgetItem, QPushButton, QVBoxLayout,
                               QWidget)

from . import __version__, theme

ASSETS = Path(__file__).parent / "assets"


class StartPage(QWidget):
    newProjectRequested = Signal()
    openProjectRequested = Signal()
    recentProjectRequested = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        outer = QVBoxLayout(self)
        outer.addStretch(2)

        logo = QSvgWidget(str(ASSETS / "logo-dark.svg"))
        logo.setFixedSize(620, 186)
        logo_row = QHBoxLayout()
        logo_row.addStretch(1)
        logo_row.addWidget(logo)
        logo_row.addStretch(1)
        outer.addLayout(logo_row)

        outer.addSpacing(28)

        buttons = QHBoxLayout()
        buttons.addStretch(1)
        new_btn = QPushButton("New project…")
        new_btn.setProperty("primary", True)
        new_btn.setFixedSize(200, 48)
        new_btn.clicked.connect(self.newProjectRequested)
        open_btn = QPushButton("Open project…")
        open_btn.setFixedSize(200, 48)
        open_btn.clicked.connect(self.openProjectRequested)
        buttons.addWidget(new_btn)
        buttons.addSpacing(16)
        buttons.addWidget(open_btn)
        buttons.addStretch(1)
        outer.addLayout(buttons)

        outer.addSpacing(32)

        recent_box = QVBoxLayout()
        caption = QLabel("Recent projects")
        caption.setProperty("dim", True)
        caption.setAlignment(Qt.AlignLeft)
        recent_box.addWidget(caption)
        self.recent_list = QListWidget()
        self.recent_list.setFixedWidth(560)
        self.recent_list.setMaximumHeight(180)
        self.recent_list.itemActivated.connect(
            lambda item: self.recentProjectRequested.emit(item.data(Qt.UserRole)))
        self.recent_list.itemClicked.connect(
            lambda item: self.recentProjectRequested.emit(item.data(Qt.UserRole)))
        recent_box.addWidget(self.recent_list)
        self.recent_empty = QLabel("No recent projects yet")
        self.recent_empty.setProperty("dim", True)
        recent_box.addWidget(self.recent_empty)

        recent_row = QHBoxLayout()
        recent_row.addStretch(1)
        recent_row.addLayout(recent_box)
        recent_row.addStretch(1)
        outer.addLayout(recent_row)

        outer.addStretch(3)

        version = QLabel(f"FSMP-kMC Studio {__version__}")
        version.setAlignment(Qt.AlignCenter)
        version.setStyleSheet(
            f"color: {theme.TEXT_DIM}; font-size: 9pt; background: transparent;")
        outer.addWidget(version)

    def set_recent(self, paths: list[str]) -> None:
        self.recent_list.clear()
        shown = 0
        for path in paths:
            p = Path(path)
            if not (p / "project.json").is_file():
                continue
            item = QListWidgetItem(f"{p.name}    —    {p}")
            item.setData(Qt.UserRole, str(p))
            self.recent_list.addItem(item)
            shown += 1
        self.recent_list.setVisible(shown > 0)
        self.recent_empty.setVisible(shown == 0)
