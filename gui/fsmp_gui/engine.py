"""Drive the engine's unit cell optimization from the GUI.

The engine already does the work: with structure = calculate and
optimize_only = true it optimizes the cell, appends an xyz frame to the
animation file after every accepted move and prints the final cell
parameters. The GUI writes a parameter file into <project>/optimize, starts
the engine as a child process there, tails the animation to show the cell
evolving live and reads the optimized cell back from the "Final params"
report.

The pure helpers (parameter text, chained-polar conversion, frame decoding)
are kept free of Qt so they are testable without a GUI session.
"""

import math
import os
import re
import shutil
from pathlib import Path

from PySide6.QtCore import QObject, QProcess, QTimer, Signal

from .molecule import Molecule
from .project import Project
from .sitemodel import SiteModel

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = "optimize"
ANIMATION = "cell_animation.xyz"
PARAMETERS = "optimize.txt"


class EngineError(Exception):
    pass


def find_engine() -> list[str] | None:
    """Command prefix that runs the engine, native first: the FSMP_ENGINE
    environment variable, a build in the repository root (fsmp.exe from a
    release or MinGW, fsmp.out from make), then fsmp on PATH. On Windows a
    Linux build is still reachable through WSL, as a last resort."""
    override = os.environ.get("FSMP_ENGINE", "")
    if override and Path(override).is_file():
        return [override]
    exe = REPO / "fsmp.exe"
    if exe.is_file():
        return [str(exe)]
    out = REPO / "fsmp.out"
    if out.is_file() and os.name != "nt":
        return [str(out)]
    found = shutil.which("fsmp")
    if found:
        return [found]
    if os.name == "nt" and out.is_file() and shutil.which("wsl"):
        drive, tail = os.path.splitdrive(str(out))
        return ["wsl", "/mnt/" + drive[0].lower() + tail.replace("\\", "/")]
    return None


# -- the unit_cell parameter (chained polar form) ---------------------------

def placements_to_chain(cell_x: float, cell_y: float, placements) -> str:
    """The unit_cell value: molecule count, the cell sides, then (r, theta,
    phi) for every molecule, each position chained from the previous one
    (the first from the cell origin)."""
    parts = [str(len(placements)), f"{cell_x:.6f}", f"{cell_y:.6f}"]
    px = py = 0.0
    for x, y, phi in placements:
        r = math.hypot(x - px, y - py)
        theta = math.degrees(math.atan2(y - py, x - px)) if r > 1e-9 else 0.0
        parts += [f"{r:.6f}", f"{theta:.6f}", f"{phi:.6f}"]
        px, py = x, y
    return "  ".join(parts)


def chain_to_placements(values) -> tuple[float, float, list]:
    """Inverse of placements_to_chain: (cell_x, cell_y, [(x, y, phi), ...])
    with the positions wrapped into the cell."""
    n = int(round(values[0]))
    cell_x, cell_y = values[1], values[2]
    out = []
    x = y = 0.0
    for i in range(n):
        r, theta, phi = values[3 + 3 * i: 6 + 3 * i]
        x += r * math.cos(math.radians(theta))
        y += r * math.sin(math.radians(theta))
        out.append((x % cell_x, y % cell_y, phi % 360.0))
    return cell_x, cell_y, out


# -- the parameter file ------------------------------------------------------

def parameter_text(potential: str, unit_cell: str) -> str:
    """Parameter file for an optimize-only run. The Monte Carlo keys are
    required by the parser but never reached: optimize_only stops the
    program right after the optimizer."""
    return (
        "# Written by FSMP-kMC Studio: unit cell optimization\n"
        f"potential = {potential}\n"
        "structure = calculate\n"
        "optimize_only = true\n"
        f"unit_cell = {unit_cell}\n"
        "molecule_model = model.xyz\n"
        f"unit_cell_name = {ANIMATION}\n"
        "sigma_mode = min_dist\n"
        "temp_from = 300\n"
        "temp_to = 300\n"
        "temp_step = 10\n"
        "um_from = 0.0\n"
        "um_to = 0.0\n"
        "um_step = 5000.0\n"
        "temperature_in_transition_zone = 900\n"
        "lambdam = 0.0\n"
        "nSteps = 1\n"
        "nStepsEq = 1\n"
        "constant_pressure = false\n"
        "kMC = true\n"
        "uc_in_x = 2\n"
        "uc_in_y = 2\n"
        "free_space = 0.24\n"
        "delta = 2.0\n"
        "delta_angle = 60.0\n"
        "widom_test_index = false\n"
        "xyz_name = trajectory.xyz\n"
    )


# -- the xyz animation --------------------------------------------------------

_LATTICE = re.compile(r'Lattice="([^"]+)"')


def decode_frame(lines: list, offsets: list, n_mols: int):
    """One xyz frame (count line, lattice line, atom rows) back into
    (cell_x, cell_y, [(x, y, phi), ...]).

    The engine writes the 3x3 tiling of the unit cell, each molecule as the
    model atoms rotated by phi and translated to the molecule position; the
    first n_mols molecules are the copy in the first tile. The center and
    the rotation are recovered from two atoms of the known model."""
    m = _LATTICE.search(lines[1])
    if m is None:
        raise ValueError("no Lattice in the frame comment")
    lat = [float(t) for t in m.group(1).split()]
    cell_x, cell_y = lat[0] / 3.0, lat[4] / 3.0
    a = len(offsets)
    # the reference pair: atom 0 and the model atom farthest from it
    j = max(range(a), key=lambda k: math.hypot(offsets[k][0] - offsets[0][0],
                                               offsets[k][1] - offsets[0][1]))
    ux, uy = offsets[j][0] - offsets[0][0], offsets[j][1] - offsets[0][1]
    out = []
    for mol in range(n_mols):
        rows = [r.split() for r in lines[2 + mol * a: 2 + (mol + 1) * a]]
        p = [(float(c[1]), float(c[2])) for c in rows]
        if math.hypot(ux, uy) > 1e-6:
            vx, vy = p[j][0] - p[0][0], p[j][1] - p[0][1]
            phi = math.degrees(math.atan2(ux * vy - uy * vx, ux * vx + uy * vy))
        else:
            phi = 0.0   # a single-site model carries no orientation
        rad = math.radians(phi)
        x = p[0][0] - (offsets[0][0] * math.cos(rad) - offsets[0][1] * math.sin(rad))
        y = p[0][1] - (offsets[0][0] * math.sin(rad) + offsets[0][1] * math.cos(rad))
        out.append((x % cell_x, y % cell_y, phi % 360.0))
    return cell_x, cell_y, out


class FrameReader:
    """Incremental reader of the growing animation file. poll() returns the
    newest complete frame decoded to (cell_x, cell_y, placements), or None
    when nothing new has appeared."""

    def __init__(self, path: Path, offsets: list, n_mols: int):
        self.path = Path(path)
        self.offsets = offsets
        self.n_mols = n_mols
        self._pos = 0
        self._tail = ""

    def poll(self):
        try:
            with open(self.path, "rb") as f:
                f.seek(self._pos)
                chunk = f.read()
                self._pos = f.tell()
        except OSError:
            return None
        if chunk:
            self._tail += chunk.decode("ascii", "replace")
        lines = self._tail.split("\n")
        newest = None
        i = 0
        while i < len(lines) - 1:   # the last element is an incomplete line
            try:
                count = int(lines[i].split()[0])
            except (ValueError, IndexError):
                return None         # not an xyz frame boundary; give up
            if i + 2 + count > len(lines) - 1:
                break               # the frame is still being written
            newest = lines[i:i + 2 + count]
            i += 2 + count
        self._tail = "\n".join(lines[i:])
        if newest is None:
            return None
        try:
            return decode_frame(newest, self.offsets, self.n_mols)
        except (ValueError, IndexError):
            return None


# -- the final report ---------------------------------------------------------

_NUMBER = re.compile(r"^Number (\d+): (\S+)", re.M)
_ENERGY = re.compile(r"^Final energy per molecule: (\S+) kJ/mol", re.M)


def final_params(text: str) -> list | None:
    """The optimized cell from the engine output: the values of the
    "Final params" listing, in order, or None."""
    parts = text.rsplit("Final params:", 1)
    if len(parts) < 2:
        return None
    values = {int(m.group(1)): float(m.group(2))
              for m in _NUMBER.finditer(parts[1])}
    if not values or sorted(values) != list(range(len(values))):
        return None
    return [values[i] for i in range(len(values))]


def final_energy(text: str) -> float | None:
    m = _ENERGY.search(text)
    return float(m.group(1)) if m else None


# -- setting up and running ----------------------------------------------------

def write_model(project: Project, run_dir: Path) -> list:
    """Write the visualization model for a run (the engine requires one)
    and return the atom offsets used to decode the animation frames. The
    atomistic model is preferred; a site model is a valid stand-in since
    only the positions matter."""
    entry = project.atomistic
    if entry is not None:
        mol = Molecule.load_xyz(project.model_path(entry))
        rows = [(a.element, a.x, a.y, a.z) for a in mol.atoms]
    elif project.site is not None:
        sm = SiteModel.load(project.model_path(project.site))
        rows = [(s.label, s.x, s.y, s.z) for s in sm.sites]
    else:
        raise EngineError("the project has no molecule model "
                          "(Molecule model tab)")
    lines = [str(len(rows)), "written by FSMP-kMC Studio for the optimizer run"]
    lines += [f"{el:<2} {x:15.8f} {y:15.8f} {z:15.8f}" for el, x, y, z in rows]
    (run_dir / "model.xyz").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return [(x, y) for _, x, y, _ in rows]


def prepare_run(project: Project, cell_x: float, cell_y: float,
                placements: list, parent=None) -> "OptimizeRun":
    """Set up an optimize-only run inside the project folder and return the
    (not yet started) run. Raises EngineError when something is missing."""
    command = find_engine()
    if command is None:
        raise EngineError("engine not found: put a release fsmp.exe in the "
                          "repository root, build one (make), or point "
                          "FSMP_ENGINE at a binary")
    if project.potential is None:
        raise EngineError("attach a potential to the project first "
                          "(Potentials tab)")
    potential = project.potential_path()
    if not potential.is_file():
        raise EngineError(f"potential file not found: {potential}")

    run_dir = project.root / RUN_DIR
    run_dir.mkdir(exist_ok=True)
    offsets = write_model(project, run_dir)
    try:
        rel = os.path.relpath(potential, run_dir).replace("\\", "/")
    except ValueError:
        raise EngineError("the potential must be on the same drive "
                          "as the project")
    chain = placements_to_chain(cell_x, cell_y, placements)
    (run_dir / PARAMETERS).write_text(parameter_text(rel, chain),
                                      encoding="utf-8")
    # the engine never overwrites an existing animation (it would shift the
    # name), so clear the previous run
    animation = run_dir / ANIMATION
    if animation.exists():
        animation.unlink()
    reader = FrameReader(animation, offsets, len(placements))
    return OptimizeRun(command, run_dir, reader, parent)


class OptimizeRun(QObject):
    """One optimization run as an engine child process."""

    frameReady = Signal(float, float, list)   # cell_x, cell_y, [(x, y, phi)]
    progress = Signal(str)                    # one line of engine output
    finished = Signal(bool, object, str)      # ok, result like frameReady, message

    def __init__(self, command: list, run_dir: Path, reader: FrameReader,
                 parent=None):
        super().__init__(parent)
        self._command = command
        self._reader = reader
        self._output = ""
        self._stopped = False
        self._process = QProcess(self)
        self._process.setWorkingDirectory(str(run_dir))
        self._process.setProcessChannelMode(QProcess.MergedChannels)
        self._process.readyReadStandardOutput.connect(self._on_output)
        self._process.finished.connect(self._on_finished)
        self._process.errorOccurred.connect(self._on_error)
        self._timer = QTimer(self)
        self._timer.setInterval(150)
        self._timer.timeout.connect(self._poll_frame)

    def start(self) -> None:
        self._process.start(self._command[0], self._command[1:] + [PARAMETERS])
        self._timer.start()

    def stop(self) -> None:
        self._stopped = True
        self._process.kill()

    def _on_output(self) -> None:
        text = bytes(self._process.readAllStandardOutput()).decode("utf-8",
                                                                   "replace")
        self._output += text
        for line in text.splitlines():
            if line.strip():
                self.progress.emit(line.strip())

    def _poll_frame(self) -> None:
        frame = self._reader.poll()
        if frame is not None:
            self.frameReady.emit(*frame)

    def _on_error(self, error) -> None:
        # a failed start never reaches _on_finished
        if error == QProcess.FailedToStart:
            self._timer.stop()
            self.finished.emit(False, None, "cannot start the engine: "
                               + " ".join(self._command))

    def _on_finished(self, code: int, status) -> None:
        self._timer.stop()
        self._poll_frame()   # pick up the frame of the last accepted move
        if self._stopped:
            self.finished.emit(False, None, "stopped")
            return
        if status != QProcess.NormalExit or code != 0:
            tail = "\n".join(self._output.strip().splitlines()[-6:])
            self.finished.emit(False, None,
                               tail or f"engine exited with code {code}")
            return
        params = final_params(self._output)
        if params is None:
            self.finished.emit(False, None,
                               "no 'Final params' in the engine output")
            return
        energy = final_energy(self._output)
        message = (f"E = {energy} kJ/mol per molecule"
                   if energy is not None else "done")
        self.finished.emit(True, chain_to_placements(params), message)
