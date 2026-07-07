"""Placeholder page for tabs that are not implemented yet."""

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QLabel, QVBoxLayout, QWidget

from .. import theme


class PlaceholderTab(QWidget):
    def __init__(self, title: str, description: str, planned: list[str],
                 parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.addStretch(2)

        heading = QLabel(title)
        heading.setAlignment(Qt.AlignCenter)
        heading.setStyleSheet(
            f"font-size: 18pt; font-weight: 600; color: {theme.TEXT_DIM};"
            "background: transparent;")
        layout.addWidget(heading)

        text = QLabel(description)
        text.setAlignment(Qt.AlignCenter)
        text.setWordWrap(True)
        text.setStyleSheet(
            f"color: {theme.TEXT_DIM}; font-size: 11pt; background: transparent;")
        layout.addWidget(text)

        layout.addSpacing(18)
        for item in planned:
            line = QLabel(f"•  {item}")
            line.setAlignment(Qt.AlignCenter)
            line.setStyleSheet(
                f"color: {theme.TEXT_DIM}; background: transparent;")
            layout.addWidget(line)

        badge = QLabel("planned")
        badge.setAlignment(Qt.AlignCenter)
        badge.setStyleSheet(
            f"color: {theme.GOLD}; font-weight: 600; letter-spacing: 2px;"
            "background: transparent; margin-top: 24px;")
        layout.addWidget(badge)
        layout.addStretch(3)
