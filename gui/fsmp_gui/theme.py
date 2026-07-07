"""Dark theme shared by the whole application.

The palette follows the project logo: teal accent, gold highlights,
deep teal-black backgrounds.
"""

from pathlib import Path

from PySide6.QtGui import QColor, QPalette
from PySide6.QtWidgets import QApplication

ARROW_DOWN = (Path(__file__).parent / "assets" / "arrow-down.svg").as_posix()

BG_WINDOW = "#0e181a"
BG_PANEL = "#132123"
BG_RAISED = "#1a2c2e"
BG_HOVER = "#22383a"
BORDER = "#284245"
TEXT = "#e4efec"
TEXT_DIM = "#8fa8a4"
ACCENT = "#1db8ad"
ACCENT_HOVER = "#27cfc3"
ACCENT_PRESSED = "#179c93"
ACCENT_BG = "#153a3b"
GOLD = "#e0c15a"
DANGER = "#e06c75"

STYLESHEET = f"""
QWidget {{
    background: {BG_WINDOW};
    color: {TEXT};
    font-family: "Segoe UI";
    font-size: 10.5pt;
}}
QMainWindow, QDialog {{ background: {BG_WINDOW}; }}

QLabel {{ background: transparent; }}
QLabel[dim="true"] {{ color: {TEXT_DIM}; }}
QLabel[heading="true"] {{ font-size: 14pt; font-weight: 600; }}
QLabel[gold="true"] {{ color: {GOLD}; font-weight: 600; }}

QPushButton {{
    background: {BG_RAISED};
    border: 1px solid {BORDER};
    border-radius: 6px;
    padding: 6px 16px;
}}
QPushButton:hover {{ background: {BG_HOVER}; border-color: {ACCENT}; }}
QPushButton:pressed {{ background: {BG_PANEL}; }}
QPushButton:disabled {{ color: {TEXT_DIM}; background: {BG_PANEL}; border-color: {BORDER}; }}
QPushButton[primary="true"] {{
    background: {ACCENT};
    border: none;
    color: #06201e;
    font-weight: 600;
}}
QPushButton[primary="true"]:hover {{ background: {ACCENT_HOVER}; }}
QPushButton[primary="true"]:pressed {{ background: {ACCENT_PRESSED}; }}
QPushButton[primary="true"]:disabled {{ background: {BG_RAISED}; color: {TEXT_DIM}; }}

QToolButton {{
    background: transparent;
    border: 1px solid transparent;
    border-radius: 6px;
    padding: 5px 10px;
}}
QToolButton:hover {{ background: {BG_HOVER}; }}
QToolButton:checked {{
    background: {ACCENT_BG};
    border-color: {ACCENT};
    color: {ACCENT_HOVER};
}}

QTabWidget::pane {{
    border: 1px solid {BORDER};
    border-radius: 8px;
    top: -1px;
}}
QTabBar::tab {{
    background: transparent;
    color: {TEXT_DIM};
    padding: 9px 22px;
    border: none;
    border-bottom: 2px solid transparent;
    margin-right: 2px;
}}
QTabBar::tab:hover {{ color: {TEXT}; }}
QTabBar::tab:selected {{
    color: {ACCENT_HOVER};
    border-bottom: 2px solid {ACCENT};
    font-weight: 600;
}}

QLineEdit, QComboBox, QDoubleSpinBox, QSpinBox {{
    background: {BG_RAISED};
    border: 1px solid {BORDER};
    border-radius: 6px;
    padding: 5px 8px;
    selection-background-color: {ACCENT};
    selection-color: #06201e;
}}
QLineEdit:focus, QComboBox:focus, QDoubleSpinBox:focus, QSpinBox:focus {{
    border-color: {ACCENT};
}}
QComboBox::drop-down {{ border: none; width: 22px; }}
QComboBox::down-arrow {{
    image: url("{ARROW_DOWN}");
    width: 10px;
    height: 6px;
    margin-right: 6px;
}}
QComboBox QAbstractItemView {{
    background: {BG_RAISED};
    border: 1px solid {BORDER};
    selection-background-color: {ACCENT_BG};
    selection-color: {TEXT};
}}

QTableView {{
    background: {BG_PANEL};
    alternate-background-color: {BG_RAISED};
    border: 1px solid {BORDER};
    border-radius: 6px;
    gridline-color: {BORDER};
    selection-background-color: {ACCENT_BG};
    selection-color: {TEXT};
}}
QHeaderView::section {{
    background: {BG_RAISED};
    color: {TEXT_DIM};
    border: none;
    border-bottom: 1px solid {BORDER};
    border-right: 1px solid {BORDER};
    padding: 5px 8px;
    font-weight: 600;
}}
QTableCornerButton::section {{ background: {BG_RAISED}; border: none; }}

QListWidget {{
    background: {BG_PANEL};
    border: 1px solid {BORDER};
    border-radius: 6px;
    padding: 3px;
}}
QListWidget::item {{
    padding: 6px 8px;
    border-radius: 4px;
}}
QListWidget::item:hover {{ background: {BG_HOVER}; }}
QListWidget::item:selected {{ background: {ACCENT_BG}; color: {TEXT}; }}

QMenuBar {{ background: {BG_WINDOW}; border-bottom: 1px solid {BORDER}; }}
QMenuBar::item {{ padding: 6px 12px; background: transparent; }}
QMenuBar::item:selected {{ background: {BG_HOVER}; border-radius: 4px; }}
QMenu {{
    background: {BG_RAISED};
    border: 1px solid {BORDER};
    border-radius: 6px;
    padding: 4px;
}}
QMenu::item {{ padding: 6px 24px 6px 16px; border-radius: 4px; }}
QMenu::item:selected {{ background: {ACCENT_BG}; }}
QMenu::separator {{ height: 1px; background: {BORDER}; margin: 4px 8px; }}

QStatusBar {{ background: {BG_PANEL}; border-top: 1px solid {BORDER}; color: {TEXT_DIM}; }}
QStatusBar::item {{ border: none; }}

QSplitter::handle {{ background: {BORDER}; }}
QSplitter::handle:horizontal {{ width: 1px; }}
QSplitter::handle:vertical {{ height: 1px; }}

QScrollBar:vertical {{
    background: transparent; width: 10px; margin: 2px;
}}
QScrollBar::handle:vertical {{
    background: {BORDER}; border-radius: 4px; min-height: 30px;
}}
QScrollBar::handle:vertical:hover {{ background: {ACCENT}; }}
QScrollBar:horizontal {{
    background: transparent; height: 10px; margin: 2px;
}}
QScrollBar::handle:horizontal {{
    background: {BORDER}; border-radius: 4px; min-width: 30px;
}}
QScrollBar::handle:horizontal:hover {{ background: {ACCENT}; }}
QScrollBar::add-line, QScrollBar::sub-line {{ height: 0; width: 0; }}
QScrollBar::add-page, QScrollBar::sub-page {{ background: transparent; }}

QToolTip {{
    background: {BG_RAISED};
    color: {TEXT};
    border: 1px solid {BORDER};
    padding: 4px 8px;
}}

QGroupBox {{
    border: 1px solid {BORDER};
    border-radius: 8px;
    margin-top: 12px;
    padding-top: 10px;
    background: transparent;
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 10px;
    padding: 0 4px;
    color: {TEXT_DIM};
    font-weight: 600;
}}

QPlainTextEdit {{
    background: {BG_PANEL};
    border: 1px solid {BORDER};
    border-radius: 6px;
    color: {TEXT};
    font-family: "Cascadia Mono", "Consolas";
    font-size: 9.5pt;
    selection-background-color: {ACCENT};
    selection-color: #06201e;
}}

QFrame[panel="true"] {{
    background: {BG_PANEL};
    border: 1px solid {BORDER};
    border-radius: 8px;
}}
"""


def apply_theme(app: QApplication) -> None:
    app.setStyle("Fusion")
    pal = QPalette()
    pal.setColor(QPalette.Window, QColor(BG_WINDOW))
    pal.setColor(QPalette.WindowText, QColor(TEXT))
    pal.setColor(QPalette.Base, QColor(BG_PANEL))
    pal.setColor(QPalette.AlternateBase, QColor(BG_RAISED))
    pal.setColor(QPalette.Text, QColor(TEXT))
    pal.setColor(QPalette.Button, QColor(BG_RAISED))
    pal.setColor(QPalette.ButtonText, QColor(TEXT))
    pal.setColor(QPalette.Highlight, QColor(ACCENT))
    pal.setColor(QPalette.HighlightedText, QColor("#06201e"))
    pal.setColor(QPalette.ToolTipBase, QColor(BG_RAISED))
    pal.setColor(QPalette.ToolTipText, QColor(TEXT))
    pal.setColor(QPalette.PlaceholderText, QColor(TEXT_DIM))
    pal.setColor(QPalette.Disabled, QPalette.Text, QColor(TEXT_DIM))
    pal.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(TEXT_DIM))
    app.setPalette(pal)
    app.setStyleSheet(STYLESHEET)
