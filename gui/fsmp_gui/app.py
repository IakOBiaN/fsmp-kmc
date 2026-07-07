"""Application entry point:  python -m fsmp_gui  or the fsmp-gui script."""

import sys

from PySide6.QtCore import QCoreApplication
from PySide6.QtWidgets import QApplication

from .main_window import MainWindow
from .theme import apply_theme


def main() -> None:
    QCoreApplication.setOrganizationName("fsmp-kmc")
    QCoreApplication.setApplicationName("fsmp-gui")
    app = QApplication(sys.argv)
    apply_theme(app)
    window = MainWindow()
    window.resize(1360, 860)
    window.show()
    sys.exit(app.exec())
