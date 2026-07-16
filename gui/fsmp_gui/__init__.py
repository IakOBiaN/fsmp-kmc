"""FSMP-kMC Studio: desktop workbench for the FSMP-kMC engine."""

import re
import sys
from pathlib import Path


def _project_version() -> str:
    """version.h is the single source of the project version: the
    repository root when running from source, baked into _internal by
    studio.spec in a frozen build."""
    if getattr(sys, "frozen", False):
        root = Path(getattr(sys, "_MEIPASS", "."))
    else:
        root = Path(__file__).resolve().parents[2]
    try:
        match = re.search(r'#define\s+FSMP_VERSION\s+"([^"]+)"',
                          (root / "version.h").read_text(encoding="utf-8"))
    except OSError:
        match = None
    return match.group(1) if match else "unknown"


__version__ = _project_version()
