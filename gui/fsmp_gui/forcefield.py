"""Reader for the v2 binary forcefield header.

The layout mirrors tools/pack_forcefield.cpp: 64 bytes, little-endian:
magic "FSMP", u32 version, u32 dtype (0 double / 1 float), u32 n_dist,
u32 n_ang, f64 min_dist, f64 dr, f64 da, f64 fold (360 = none), zero pad.
"""

import struct
from dataclasses import dataclass
from pathlib import Path

HEADER_BYTES = 64
MAGIC = b"FSMP"


class ForcefieldError(Exception):
    pass


def _human_size(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{n} B"
        n /= 1024.0
    return f"{n:.1f} GB"


@dataclass
class ForcefieldInfo:
    path: Path
    version: int
    dtype: str  # "double" or "float"
    n_dist: int
    n_ang: int
    min_dist: float
    dr: float
    da: float
    fold: float  # 360 = no folding
    file_size: int

    @property
    def r_max(self) -> float:
        return self.min_dist + (self.n_dist - 1) * self.dr

    @property
    def n_values(self) -> int:
        return self.n_dist * self.n_ang * self.n_ang

    @property
    def data_bytes(self) -> int:
        return self.n_values * (4 if self.dtype == "float" else 8)

    def summary(self) -> str:
        fold = "none" if self.fold >= 359.999 else f"{self.fold:g}°"
        return (f"{self.dtype}, {self.n_dist} × {self.n_ang} × {self.n_ang} points, "
                f"fold {fold}\n"
                f"r = {self.min_dist:g} … {self.r_max:g} Å (dr {self.dr:g} Å, "
                f"da {self.da:g}°), {_human_size(self.file_size)}")


def read_header(path: str | Path) -> ForcefieldInfo:
    path = Path(path)
    try:
        with path.open("rb") as f:
            data = f.read(HEADER_BYTES)
        file_size = path.stat().st_size
    except OSError as e:
        raise ForcefieldError(str(e))
    if len(data) < HEADER_BYTES:
        raise ForcefieldError("file is shorter than the 64-byte header")
    if data[:4] != MAGIC:
        raise ForcefieldError("not an FSMP forcefield (bad magic)")
    version, dtype, n_dist, n_ang = struct.unpack_from("<4I", data, 4)
    min_dist, dr, da, fold = struct.unpack_from("<4d", data, 20)
    if version != 2:
        raise ForcefieldError(f"unsupported format version {version} (expected 2)")
    if dtype not in (0, 1):
        raise ForcefieldError(f"unknown dtype {dtype}")
    info = ForcefieldInfo(path, version, "float" if dtype else "double",
                          n_dist, n_ang, min_dist, dr, da, fold, file_size)
    expected = HEADER_BYTES + info.data_bytes
    if file_size != expected:
        raise ForcefieldError(f"file size mismatch: {file_size} bytes on disk, "
                              f"the header implies {expected}")
    return info
