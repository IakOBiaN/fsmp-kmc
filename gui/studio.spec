# PyInstaller build of FSMP-kMC Studio: one folder (dist/studio) with a
# windowed executable. Build from the repository root:
#
#   pyinstaller --noconfirm --distpath gui/dist --workpath gui/build gui/studio.spec
#
# The release workflow copies the folder's contents next to the engine
# binaries and the configs/ and samples/ folders; the app finds them there
# through fsmp_gui.engine.app_root().

import re
import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_all

# RDKit ships compiled libraries and data files (atom typing, MMFF tables)
# that PyInstaller does not pick up automatically
rd_datas, rd_binaries, rd_hidden = collect_all("rdkit")

# Version resource for the Windows executable (Explorer Properties ->
# Details), from version.h, the single project-version source
VERSION = re.search(r'#define\s+FSMP_VERSION\s+"([^"]+)"',
                    (Path(SPECPATH) / ".." / "version.h").read_text()).group(1)
version_resource = None
if sys.platform == "win32":
    from PyInstaller.utils.win32.versioninfo import (
        FixedFileInfo, StringFileInfo, StringStruct, StringTable,
        VarFileInfo, VarStruct, VSVersionInfo)
    quad = tuple(int(n) for n in (VERSION.split(".") + ["0"] * 4)[:4])
    version_resource = VSVersionInfo(
        ffi=FixedFileInfo(filevers=quad, prodvers=quad),
        kids=[StringFileInfo([StringTable("040904B0", [
                  StringStruct("CompanyName", "The FSMP-kMC authors"),
                  StringStruct("FileDescription", "FSMP-kMC Studio"),
                  StringStruct("FileVersion", VERSION),
                  StringStruct("InternalName", "FSMP-kMC Studio"),
                  StringStruct("OriginalFilename", "FSMP-kMC Studio.exe"),
                  StringStruct("ProductName", "FSMP-kMC"),
                  StringStruct("ProductVersion", VERSION),
                  StringStruct("LegalCopyright",
                               "(C) The FSMP-kMC authors. GPL-3.0")])]),
              VarFileInfo([VarStruct("Translation", [0x409, 1200])])])

a = Analysis(
    ["studio_launcher.py"],
    pathex=[],
    binaries=rd_binaries,
    # version.h rides along so the frozen app knows the project version
    datas=[("fsmp_gui/assets", "fsmp_gui/assets"),
           ("../version.h", ".")] + rd_datas,
    hiddenimports=rd_hidden,
    excludes=["tkinter"],
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz,
    a.scripts,
    exclude_binaries=True,
    name="FSMP-kMC Studio",
    console=False,
    icon="studio.ico" if sys.platform == "win32" else None,
    version=version_resource,
)
coll = COLLECT(exe, a.binaries, a.datas, name="studio")

if sys.platform == "darwin":
    app = BUNDLE(coll, name="FSMP-kMC Studio.app",
                 bundle_identifier="io.github.iakobian.fsmp-kmc")
