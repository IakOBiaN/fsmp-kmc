"""Run management: the (T, u_m) point grid must mirror the engine loops,
the log watcher must track progress incrementally, and the round trip
parameter file -> engine -> statistics must work end to end (the last test
runs the real engine detached and is skipped when it is absent).

Run from the repository root:
    gui/.venv/Scripts/python gui/tests/test_runs.py
"""

import sys
import tempfile
import time
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsmp_gui.engine import find_engine
from fsmp_gui.runs import (DONE, INTERRUPTED, RUNNING, STOPPED, LogWatch,
                           create_run, is_alive, loop_points, loop_values,
                           point_at, read_params, read_statistics,
                           run_parameters, stop, waiter_main)

REPO = Path(__file__).resolve().parents[2]
GRID = REPO / "tests" / "data" / "TMA_simple_2020_s4.v2.bin"
SITE = REPO / "models" / "TMA_simplified_2020.site"
MODEL = REPO / "models" / "trimesic_acid.xyz"

FORM = {"temp_from": 300.0, "temp_to": 300.0, "temp_step": 10.0,
        "um_from": 0.0, "um_to": 0.0, "um_step": 5000.0,
        "nSteps": 3, "nStepsEq": 1, "kMC": True,
        "delta": 2.0, "delta_angle": 60.0, "widom": False,
        "sigma_mode": "manual", "sigma": 11.052,
        "constant_pressure": False, "constant_pressure_value": 0.0,
        "mask": False, "mask_free_radius": 3.0, "mask_ramp_width": 2.0,
        "mask_penalty": 25000.0, "seed": 12345}

SIM = {"uc_in_x": 2, "uc_in_y": 2, "free_space": 0.24, "lambdam": 0.0,
       "temperature_in_transition_zone": 900.0}


class TestLoops(unittest.TestCase):
    def test_ascending(self):
        self.assertEqual(loop_values(300, 350, 10), [300, 310, 320, 330,
                                                     340, 350])

    def test_descending_negative_step(self):
        # the engine takes |step| and walks down when from > to
        self.assertEqual(loop_values(350, 330, -10), [350, 340, 330])

    def test_single_point(self):
        self.assertEqual(loop_values(300, 300, 10), [300])
        self.assertEqual(loop_values(300, 305, 10), [300])
        self.assertEqual(loop_values(300, 400, 0), [300])

    def test_point_order_um_outer(self):
        ts, ums = loop_points({"temp_from": "300", "temp_to": "310",
                               "temp_step": "10", "um_from": "0",
                               "um_to": "-5000", "um_step": "5000"})
        self.assertEqual(ts, [300, 310])
        self.assertEqual(ums, [0, -5000])
        order = [point_at(ts, ums, i) for i in range(4)]
        self.assertEqual(order, [(300, 0), (310, 0), (300, -5000),
                                 (310, -5000)])


class TestLogWatch(unittest.TestCase):
    def test_incremental_progress(self):
        with tempfile.TemporaryDirectory() as td:
            log = Path(td) / "run.log"
            watch = LogWatch(log)
            self.assertFalse(watch.poll())

            log.write_text("Parameters read from run.txt\n"
                           "Unit cell optimization: 9 degrees of freedom\n"
                           "Cell scaling: factor 1, energy -60 kJ/mol\n")
            watch.poll()
            self.assertEqual(watch.points_started, 0)
            self.assertEqual(watch.state(), RUNNING)

            with open(log, "a") as f:
                f.write("_________INITIAL DATA_________\n\nu_m: 0\n"
                        "5 %\n17 %\n")
            watch.poll()
            self.assertEqual((watch.points_started, watch.percent), (1, 17))

            with open(log, "a") as f:  # second point resets the percent
                f.write("100 %\n_________INITIAL DATA_________\n7 %")
            watch.poll()
            self.assertEqual((watch.points_started, watch.percent), (2, 0))
            with open(log, "a") as f:
                f.write("\nFSMP_EXIT:0\n")
            watch.poll()
            self.assertEqual(watch.percent, 7)
            self.assertEqual(watch.state(), DONE)


class TestStatistics(unittest.TestCase):
    def test_reads_rows_after_header(self):
        text = ("Number of particles: 8\nTotal number of MCS: 3  MCS...\n"
                "Maximal displacement, A: 2\nLambda0: 0.57 Lambdam: 0\n"
                "////////////\n\n"
                "T, K\tu_m, kJ/mol\tDensity, mkmol/m2\n"
                "300\t0\t1.55\n"
                "310\t0\tnan\n")
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "statistics.dat"
            path.write_text(text)
            names, rows = read_statistics(path)
        self.assertEqual(names, ["T, K", "u_m, kJ/mol", "Density, mkmol/m2"])
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0], [300.0, 0.0, 1.55])
        self.assertNotEqual(rows[1][2], rows[1][2])   # nan


class TestParameters(unittest.TestCase):
    def test_round_trip_and_required_keys(self):
        text = run_parameters(FORM, "../grid.bin", "1 10 10 0 0 84", SIM)
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "run.txt"
            path.write_text(text)
            params = read_params(path)
        for key in ("potential", "structure", "sigma_mode", "sigma",
                    "temp_from", "temp_to", "temp_step", "um_from", "um_to",
                    "um_step", "temperature_in_transition_zone", "lambdam",
                    "nSteps", "nStepsEq", "constant_pressure", "kMC",
                    "uc_in_x", "uc_in_y", "free_space", "molecule_model",
                    "delta", "delta_angle", "widom_test_index",
                    "unit_cell_name", "xyz_name", "unit_cell", "seed"):
            self.assertIn(key, params, msg=key)
        self.assertEqual(params["structure"], "calculate")
        self.assertEqual(params["kMC"], "true")
        self.assertNotIn("stabilization_mask", params)
        self.assertNotIn("constant_pressure_value", params)

    def test_optional_blocks(self):
        form = dict(FORM, mask=True, constant_pressure=True, seed=0,
                    sigma_mode="min_dist")
        params_text = run_parameters(form, "g.bin", "1 10 10 0 0 0", SIM)
        self.assertIn("stabilization_mask = true", params_text)
        self.assertIn("mask_penalty = 25000", params_text)
        self.assertIn("constant_pressure_value = 0", params_text)
        self.assertNotIn("seed", params_text)
        self.assertNotIn("sigma =", params_text)


class TestLegacyWslRun(unittest.TestCase):
    def test_state_comes_from_files_only(self):
        """A run made through WSL before the native era carries a Linux
        PID: it must never be probed or killed as a Windows PID."""
        with tempfile.TemporaryDirectory() as td:
            run_dir = Path(td)
            (run_dir / "run.json").write_text('{"label": "old", "kind": "wsl"}')
            (run_dir / "engine.pid").write_text("12345")
            (run_dir / "run.log").write_text(
                "_________INITIAL DATA_________\n5 %\n")
            self.assertFalse(is_alive(run_dir))
            stop(run_dir)   # a no-op, not a TerminateProcess on pid 12345
            watch = LogWatch(run_dir / "run.log")
            watch.poll()
            self.assertEqual(watch.state(alive=is_alive(run_dir)),
                             INTERRUPTED)


class TestWaiter(unittest.TestCase):
    def test_records_pid_log_and_exit_marker(self):
        """waiter_main with a python script standing in for the engine."""
        import os
        with tempfile.TemporaryDirectory() as td:
            script = Path(td) / "engine.py"
            script.write_text("print('engine says hi')")
            cwd = os.getcwd()
            os.chdir(td)
            try:
                waiter_main(sys.executable, str(script))
            finally:
                os.chdir(cwd)
            log = (Path(td) / "run.log").read_text()
            self.assertIn("engine says hi", log)
            self.assertTrue(log.rstrip().endswith("FSMP_EXIT:0"))
            self.assertGreater(int((Path(td) / "engine.pid").read_text()), 0)


def _e2e_project(td):
    from fsmp_gui.molecule import Molecule
    from fsmp_gui.project import Project
    from fsmp_gui.sitemodel import SiteModel

    project = Project.create(Path(td) / "proj", "e2e")
    project.set_atomistic("TMA", Molecule.load_xyz(MODEL))
    project.set_site("TMA", SiteModel.load(SITE))
    project.set_potential("small-grid", GRID)
    project.set_unit_cell(11.6, 19.8, [
        {"x": 0.0, "y": 0.0, "phi": 84.0},
        {"x": 5.8, "y": 9.9, "phi": 96.0}])
    project.set_simulation_cell(SIM)
    return project


@unittest.skipUnless(find_engine() and GRID.is_file() and SITE.is_file(),
                     "engine or test grid not present")
class TestEndToEnd(unittest.TestCase):
    def test_detached_run_completes(self):

        # ignore_cleanup_errors: the detached waiter can hold the log
        # handle a moment after the exit marker appears
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
            run_dir = create_run(_e2e_project(td), "e2e-run", FORM)
            watch = LogWatch(run_dir / "run.log")
            deadline = time.time() + 300
            while time.time() < deadline:
                watch.poll()
                if watch.exit_code is not None:
                    break
                time.sleep(0.5)
            self.assertEqual(watch.exit_code, 0,
                             msg=(run_dir / "run.log").read_text()[-2000:])
            self.assertEqual(watch.state(), DONE)
            self.assertEqual(watch.points_started, 1)
            self.assertTrue((run_dir / "engine.pid").is_file())

            names, rows = read_statistics(run_dir / "statistics.dat")
            self.assertEqual(names[0], "T, K")
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0][0], 300.0)   # the requested temperature

            # the trajectory indexes and its frames parse
            from fsmp_gui.trajectory import FrameIndexer, read_frame
            got = []
            indexer = FrameIndexer(run_dir / "trajectory.xyz")
            indexer.indexed.connect(lambda off, end: got.append((off, end)))
            indexer.run()   # synchronously, without the thread machinery
            self.assertTrue(got and got[0][0], msg="no frames indexed")
            lattice, atoms = read_frame(run_dir / "trajectory.xyz",
                                        got[0][0][-1])
            self.assertIsNotNone(lattice)
            self.assertEqual(len(atoms), 2 * SIM["uc_in_x"] * SIM["uc_in_y"]
                             * 21)   # molecules x atoms of the TMA model

    def test_stop_terminates_the_engine(self):
        # a run long enough to still be going when we shoot it down
        form = dict(FORM, nSteps=2000000, nStepsEq=1000000)
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
            run_dir = create_run(_e2e_project(td), "e2e-stop", form)
            # engine.pid may exist before the pid is flushed into it
            deadline = time.time() + 60
            while time.time() < deadline and not is_alive(run_dir):
                time.sleep(0.2)
            self.assertTrue(is_alive(run_dir))

            stop(run_dir)
            watch = LogWatch(run_dir / "run.log")
            deadline = time.time() + 30
            while time.time() < deadline:
                watch.poll()
                if watch.exit_code is not None:
                    break
                time.sleep(0.5)
            self.assertIsNotNone(watch.exit_code, msg="no exit marker")
            self.assertEqual(watch.state(), STOPPED)
            self.assertFalse(is_alive(run_dir))


if __name__ == "__main__":
    unittest.main(verbosity=2)
