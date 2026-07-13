"""Unit cell file input/output (.cell): a plain xyz-like text format, so
cells can be saved, shared and loaded like molecule models. The reference
structures from StructureGenerator.h live in the repository cells/ folder
in this format.

    <n_molecules>
    <cell_x> <cell_y> [comment]
    x  y  phi        one line per molecule; angstroms and degrees
"""

from pathlib import Path

SUFFIX = ".cell"


def load_cell(path: str | Path):
    """(cell_x, cell_y, [(x, y, phi), ...], comment)"""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError("empty file")
    try:
        count = int(lines[0].split()[0])
    except (ValueError, IndexError):
        raise ValueError("line 1: expected the number of molecules")
    head = lines[1].split() if len(lines) > 1 else []
    try:
        cell_x, cell_y = float(head[0]), float(head[1])
    except (ValueError, IndexError):
        raise ValueError("line 2: expected 'cell_x cell_y [comment]'")
    comment = " ".join(head[2:])
    placements = []
    for i in range(count):
        ln = 2 + i
        if ln >= len(lines):
            raise ValueError(f"expected {count} molecules, file ends after {i}")
        cols = lines[ln].split()
        if len(cols) < 3:
            raise ValueError(f"line {ln + 1}: expected 'x y phi'")
        try:
            placements.append((float(cols[0]), float(cols[1]),
                               float(cols[2])))
        except ValueError:
            raise ValueError(f"line {ln + 1}: bad number")
    return cell_x, cell_y, placements, comment


def save_cell(path: str | Path, cell_x: float, cell_y: float,
              placements, comment: str = "") -> None:
    lines = [str(len(placements)),
             f"{cell_x:.6f} {cell_y:.6f}" + (f"  {comment}" if comment else "")]
    for x, y, phi in placements:
        lines.append(f"{x:12.6f} {y:12.6f} {phi:10.4f}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
