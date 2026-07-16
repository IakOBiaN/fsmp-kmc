"""Production runs for the Run tab. Each run lives in its own folder
<project>/runs/<name>/ with everything it needs and everything it produces:

    run.txt         the engine parameter file (the run's true metadata)
    run.json        label, launcher kind, creation time
    model.xyz       the visualization molecule model
    run.log         engine stdout+stderr; the launcher appends FSMP_EXIT:<code>
    engine.pid      PID of the engine process
    statistics.dat  one row per (T, u_m) point
    trajectory.xyz, unit_cell.xyz

The engine is launched detached from the GUI, so closing the window does not
kill multi-day runs; the GUI recovers every run's state from these files.
No Qt in this module: everything is testable standalone.
"""

import json
import os
import re
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from .engine import find_engine, placements_to_chain, write_model
from .project import Project, safe_filename

RUNS_DIR = "runs"
PARAMS = "run.txt"
META = "run.json"
LOG = "run.log"
PIDFILE = "engine.pid"
STATISTICS = "statistics.dat"
TRAJECTORY = "trajectory.xyz"
CELL_ANIMATION = "unit_cell.xyz"
EXIT_MARK = "FSMP_EXIT:"

# CREATE_NO_WINDOW alone: combined with DETACHED_PROCESS Windows ignores it
# and the launcher gets a visible console. Children survive the GUI anyway.
_NO_WINDOW = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0


class RunError(Exception):
    pass


# -- the (T, u_m) point grid --------------------------------------------------

def loop_values(v_from: float, v_to: float, v_step: float) -> list:
    """The values a program_body.cpp while-loop visits: repeated addition of
    |step| from v_from towards v_to, endpoints as the engine treats them."""
    v_step = abs(v_step)
    out = [v_from]
    if v_step <= 0:
        return out
    v = v_from
    if v_from < v_to:
        while True:
            v += v_step
            if v > v_to:
                break
            out.append(v)
    else:
        while True:
            v -= v_step
            if v < v_to:
                break
            out.append(v)
    return out


def loop_points(params: dict) -> tuple[list, list]:
    """(temperatures, u_m values) of a run; u_m is the outer loop."""
    ts = loop_values(float(params["temp_from"]), float(params["temp_to"]),
                     float(params["temp_step"]))
    ums = loop_values(float(params["um_from"]), float(params["um_to"]),
                      float(params["um_step"]))
    return ts, ums


def point_at(ts: list, ums: list, index: int) -> tuple[float, float]:
    """(T, u_m) of the 0-based point index (T is the inner loop)."""
    return ts[index % len(ts)], ums[index // len(ts)]


# -- the parameter file --------------------------------------------------------

def read_params(path: Path) -> dict:
    """key = value pairs of a parameter file (comments stripped)."""
    out = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            key, value = line.split("=", 1)
            out[key.strip()] = value.strip()
    return out


def run_parameters(form: dict, potential: str, unit_cell: str,
                   sim: dict) -> str:
    """The parameter file of a production run: the project potential and
    unit cell, the simulation-cell settings from tab 5 and the run form."""
    flag = lambda v: "true" if v else "false"
    lines = [
        "# Written by FSMP-kMC Studio: production run",
        f"potential = {potential}",
        "structure = calculate",
        f"unit_cell = {unit_cell}",
        "molecule_model = model.xyz",
        f"uc_in_x = {int(sim['uc_in_x'])}",
        f"uc_in_y = {int(sim['uc_in_y'])}",
        f"free_space = {sim['free_space']:g}",
        f"temperature_in_transition_zone = {sim['temperature_in_transition_zone']:g}",
        f"lambdam = {sim['lambdam']:g}",
        f"temp_from = {form['temp_from']:g}",
        f"temp_to = {form['temp_to']:g}",
        f"temp_step = {form['temp_step']:g}",
        f"um_from = {form['um_from']:g}",
        f"um_to = {form['um_to']:g}",
        f"um_step = {form['um_step']:g}",
        f"nSteps = {int(form['nSteps'])}",
        f"nStepsEq = {int(form['nStepsEq'])}",
        f"kMC = {flag(form['kMC'])}",
        f"delta = {form['delta']:g}",
        f"delta_angle = {form['delta_angle']:g}",
        f"widom_test_index = {flag(form['widom'])}",
        f"constant_pressure = {flag(form['constant_pressure'])}",
    ]
    if form["constant_pressure"]:
        lines.append(f"constant_pressure_value = {form['constant_pressure_value']:g}")
    if form["sigma_mode"] == "manual":
        lines += ["sigma_mode = manual", f"sigma = {form['sigma']:g}"]
    else:
        lines.append(f"sigma_mode = {form['sigma_mode']}")
    if form.get("mask"):
        lines += ["stabilization_mask = true",
                  f"mask_free_radius = {form['mask_free_radius']:g}",
                  f"mask_ramp_width = {form['mask_ramp_width']:g}",
                  f"mask_penalty = {form['mask_penalty']:g}"]
    if form.get("seed"):
        lines.append(f"seed = {int(form['seed'])}")
    lines += [f"unit_cell_name = {CELL_ANIMATION}",
              f"xyz_name = {TRAJECTORY}",
              f"statistics_name = {STATISTICS}"]
    return "\n".join(lines) + "\n"


# -- launching and stopping -----------------------------------------------------

# The detached wrappers stay referenced here: a garbage-collected Popen of a
# still-running child (which a wrapper is, by design) raises ResourceWarning.
_launched: list = []


def waiter_main(engine: str, params: str) -> None:
    """Body of the frozen app's detached launcher: studio_launcher.py calls
    this when the executable is re-entered with --run-waiter (a frozen
    sys.executable is the Studio itself, so launch() cannot hand it the
    `python -c` wrapper below). Same contract as that wrapper: the engine
    PID to engine.pid, output to run.log, the exit marker appended, a
    negative wait() folded to 128+N."""
    with open(LOG, "w") as log:
        proc = subprocess.Popen([engine, params], stdout=log,
                                stderr=subprocess.STDOUT,
                                stdin=subprocess.DEVNULL,
                                creationflags=_NO_WINDOW)
    Path(PIDFILE).write_text(str(proc.pid))
    code = proc.wait()
    if code < 0:
        code = 128 - code
    with open(LOG, "a") as log:
        log.write(f"\n{EXIT_MARK}{code}\n")


def launch(run_dir: Path, command: list) -> str:
    """Start the engine detached from the GUI. Returns the launcher kind.
    The waiter writes the engine PID to engine.pid, sends all output to
    run.log and appends FSMP_EXIT:<code> when the engine finishes."""
    if getattr(sys, "frozen", False):
        # a frozen sys.executable is the Studio itself; its launcher
        # recognizes --run-waiter and calls waiter_main
        args = [sys.executable, "--run-waiter", command[0], PARAMS]
    else:
        # The engine gets CREATE_NO_WINDOW of its own: the wrapper has no
        # console, so a console child would otherwise open a visible one.
        # A negative wait() status (POSIX signal death) is folded to the
        # 128+N shell convention, so state() sees 143 for SIGTERM.
        wrapper = (
            "import os, subprocess\n"
            "flags = 0x08000000 if os.name == 'nt' else 0\n"
            f"p = subprocess.Popen({command + [PARAMS]!r},"
            f" stdout=open({LOG!r}, 'w'), stderr=subprocess.STDOUT,"
            " creationflags=flags)\n"
            f"open({PIDFILE!r}, 'w').write(str(p.pid))\n"
            "code = p.wait()\n"
            "if code < 0:\n"
            "    code = 128 - code\n"
            f"open({LOG!r}, 'a').write('\\n{EXIT_MARK}' + str(code) + '\\n')\n")
        args = [sys.executable, "-c", wrapper]
    _launched.append(
        subprocess.Popen(args, cwd=str(run_dir), creationflags=_NO_WINDOW,
                         stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL))
    return "native"


def _pid(run_dir: Path) -> int | None:
    try:
        return int((run_dir / PIDFILE).read_text().split()[0])
    except (OSError, ValueError, IndexError):
        return None


def _native(run_dir: Path) -> bool:
    """True for runs made by this launcher. Runs launched through WSL
    before the native era carry kind "wsl" and a Linux PID: their state
    comes from the run files only, and that PID must never be checked or
    killed in the Windows pid namespace."""
    try:
        meta = json.loads((run_dir / META).read_text(encoding="utf-8"))
        return meta["kind"] == "native"
    except (OSError, ValueError, KeyError):
        return False


_PROCESS_TERMINATE = 0x0001
_PROCESS_QUERY_LIMITED = 0x1000
_STILL_ACTIVE = 259


def _win_alive(pid: int) -> bool:
    import ctypes
    kernel = ctypes.windll.kernel32
    handle = kernel.OpenProcess(_PROCESS_QUERY_LIMITED, False, pid)
    if not handle:
        return False
    try:
        code = ctypes.c_ulong()
        if not kernel.GetExitCodeProcess(handle, ctypes.byref(code)):
            return False
        return code.value == _STILL_ACTIVE
    finally:
        kernel.CloseHandle(handle)


def _win_kill(pid: int, code: int) -> None:
    """TerminateProcess with the given exit code: Windows has no graceful
    signal, but the code makes the wrapper's exit marker say why it died."""
    import ctypes
    kernel = ctypes.windll.kernel32
    handle = kernel.OpenProcess(_PROCESS_TERMINATE, False, pid)
    if handle:
        kernel.TerminateProcess(handle, code)
        kernel.CloseHandle(handle)


def stop(run_dir: Path) -> None:
    """Terminate the engine. POSIX asks with SIGTERM and escalates to
    SIGKILL; Windows terminates with exit code 143 so the run reads as
    stopped, the same as a SIGTERM death."""
    pid = _pid(run_dir)
    if pid is None or not _native(run_dir):
        return
    if os.name == "nt":
        _win_kill(pid, 143)
    else:
        try:
            os.kill(pid, signal.SIGTERM)
        except OSError:
            return
        for _ in range(6):
            time.sleep(0.5)
            if not is_alive(run_dir):
                return
        try:
            os.kill(pid, signal.SIGKILL)
        except OSError:
            pass


def is_alive(run_dir: Path) -> bool:
    pid = _pid(run_dir)
    if pid is None or not _native(run_dir):
        return False
    if os.name == "nt":
        return _win_alive(pid)
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except OSError:
        return True   # exists but owned by another user


# -- run state from the log -----------------------------------------------------

_PERCENT = re.compile(r"^(\d+) %\s*$")

RUNNING = "running"
DONE = "done"
FAILED = "failed"
STOPPED = "stopped"
INTERRUPTED = "interrupted"


class LogWatch:
    """Incremental reader of run.log keeping the progress counters:
    points_started (INITIAL DATA markers), percent within the current point
    and the exit code once the FSMP_EXIT line appears."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.points_started = 0
        self.percent = 0
        self.exit_code: int | None = None
        self._pos = 0
        self._tail = ""

    def poll(self) -> bool:
        """Consume newly appended text; True when anything changed."""
        try:
            with open(self.path, "rb") as f:
                f.seek(self._pos)
                chunk = f.read()
                self._pos = f.tell()
        except OSError:
            return False
        if not chunk:
            return False
        self._tail += chunk.decode("utf-8", "replace")
        lines = self._tail.split("\n")
        self._tail = lines.pop()          # keep the incomplete last line
        for line in lines:
            line = line.strip()
            if "_________INITIAL DATA_________" in line:
                self.points_started += 1
                self.percent = 0
            elif line.startswith(EXIT_MARK):
                try:
                    self.exit_code = int(line[len(EXIT_MARK):])
                except ValueError:
                    pass
            else:
                m = _PERCENT.match(line)
                if m:
                    self.percent = int(m.group(1))
        return True

    def state(self, alive: bool | None = None) -> str:
        """The run state; pass alive= from is_alive when the exit marker is
        absent to distinguish a running engine from an interrupted one."""
        if self.exit_code is None:
            return RUNNING if alive in (True, None) else INTERRUPTED
        if self.exit_code == 0:
            return DONE
        if self.exit_code in (130, 137, 143):   # SIGINT / SIGKILL / SIGTERM
            return STOPPED
        return FAILED


# -- statistics ------------------------------------------------------------------

def read_statistics(path: Path) -> tuple[list, list]:
    """(column names, rows) of the engine statistics file. The header block
    is skipped; the column-name line starts with "T, K"."""
    try:
        text = Path(path).read_text(encoding="utf-8")
    except OSError:
        return [], []
    names, rows = [], []
    for line in text.splitlines():
        if line.startswith("T, K"):
            names = [n.strip() for n in line.split("\t")]
        elif names and line.strip():
            try:
                row = [float(t) for t in line.split("\t")]
            except ValueError:
                continue
            if len(row) == len(names):
                rows.append(row)
    return names, rows


# -- creating and listing runs -----------------------------------------------------

def create_run(project: Project, label: str, form: dict) -> Path:
    """Prepare the run folder, write everything the engine needs and launch
    it detached. Raises RunError when a project piece is missing."""
    command = find_engine()
    if command is None:
        raise RunError("engine not found: build fsmp.out in the repository "
                       "root (make) or put fsmp.exe next to it")
    if project.potential is None:
        raise RunError("attach a potential to the project first (Potentials tab)")
    potential = project.potential_path()
    if not potential.is_file():
        raise RunError(f"potential file not found: {potential}")
    uc = project.unit_cell
    if uc is None or not uc["molecules"]:
        raise RunError("save a unit cell on the Unit cell tab first")
    sim = project.simulation_cell
    if sim is None:
        raise RunError("save the simulation cell on the Simulation cell tab first")

    root = project.root / RUNS_DIR
    root.mkdir(exist_ok=True)
    slug = safe_filename(label) or "run"
    run_dir = root / slug
    n = 2
    while run_dir.exists():
        run_dir = root / f"{slug}-{n}"
        n += 1
    run_dir.mkdir()

    write_model(project, run_dir)
    try:
        rel = os.path.relpath(potential, run_dir).replace("\\", "/")
    except ValueError:
        raise RunError("the potential must be on the same drive as the project")
    chain = placements_to_chain(uc["cell_x"], uc["cell_y"],
                                [(m["x"], m["y"], m["phi"])
                                 for m in uc["molecules"]])
    (run_dir / PARAMS).write_text(run_parameters(form, rel, chain, sim),
                                  encoding="utf-8")
    (run_dir / META).write_text(json.dumps({
        "label": label, "kind": "native",
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds")},
        indent=2), encoding="utf-8")
    launch(run_dir, command)
    return run_dir


def run_label(run_dir: Path) -> str:
    try:
        return json.loads((run_dir / META).read_text(encoding="utf-8"))["label"]
    except (OSError, ValueError, KeyError):
        return run_dir.name


def list_runs(project: Project) -> list:
    """Run folders of the project, newest first."""
    root = project.root / RUNS_DIR
    if not root.is_dir():
        return []
    dirs = [d for d in root.iterdir() if (d / PARAMS).is_file()]
    return sorted(dirs, key=lambda d: (d / PARAMS).stat().st_mtime,
                  reverse=True)
