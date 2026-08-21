import sys
from pathlib import Path

_GUI = Path(__file__).resolve().parents[1] / "gui"
if str(_GUI) not in sys.path:
    sys.path.insert(0, str(_GUI))
