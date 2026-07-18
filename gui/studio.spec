# PyInstaller build of FSMP-kMC Studio: one folder (dist/studio) with a
# windowed executable. Build from the repository root:
#
#   pyinstaller --noconfirm --distpath gui/dist --workpath gui/build gui/studio.spec
#
# The release workflow copies the folder's contents next to the engine
# binaries and the models/, cells/ and configs/ folders; the app finds them
# there through fsmp_gui.engine.app_root().

import sys

from PyInstaller.utils.hooks import collect_all

# RDKit ships compiled libraries and data files (atom typing, MMFF tables)
# that PyInstaller does not pick up automatically
rd_datas, rd_binaries, rd_hidden = collect_all("rdkit")

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
)
coll = COLLECT(exe, a.binaries, a.datas, name="studio")

if sys.platform == "darwin":
    app = BUNDLE(coll, name="FSMP-kMC Studio.app",
                 bundle_identifier="io.github.iakobian.fsmp-kmc")
