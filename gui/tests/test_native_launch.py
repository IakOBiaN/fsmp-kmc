"""The native (no WSL) detached launcher on the current platform, with a
Python stand-in for the engine: the launch/pid/log plumbing, is_alive, stop
and the exit-code mapping. The real-engine end-to-end runs live in
test_runs; these tests need no engine binary at all."""

import json
import sys
import tempfile
import time
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsmp_gui import runs
from fsmp_gui.runs import (DONE, LOG, PIDFILE, STOPPED, LogWatch, is_alive,
                           launch, stop)

SLEEPER = ("import sys, time\n"
           "print('fake engine, parameters:', sys.argv[1], flush=True)\n"
           "time.sleep(120)\n")
QUICK = "print('fake engine: finished cleanly')\n"


def wait_for(check, timeout=15.0):
    deadline = time.time() + timeout
    while time.time() < deadline:
        if check():
            return True
        time.sleep(0.1)
    return False


class NativeLaunch(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(self._tmp.cleanup)

    def _run_dir(self, script: str) -> Path:
        d = Path(self._tmp.name)
        (d / "engine.py").write_text(script, encoding="utf-8")
        (d / runs.PARAMS).write_text("# fake parameters\n", encoding="utf-8")
        (d / runs.META).write_text(
            json.dumps({"label": "t", "kind": "native"}), encoding="utf-8")
        return d

    def test_stop_kills_and_reads_as_stopped(self):
        d = self._run_dir(SLEEPER)
        kind = launch(d, [sys.executable, str(d / "engine.py")])
        self.assertEqual(kind, "native")
        self.assertTrue(wait_for(lambda: (d / PIDFILE).is_file()),
                        "no pid file")
        self.assertTrue(wait_for(lambda: is_alive(d)), "engine not alive")

        stop(d)
        self.assertTrue(wait_for(lambda: not is_alive(d)), "still alive")
        watch = LogWatch(d / LOG)

        def has_exit():
            watch.poll()
            return watch.exit_code is not None

        self.assertTrue(wait_for(has_exit), "no exit marker in the log")
        self.assertEqual(watch.exit_code, 143)
        self.assertEqual(watch.state(alive=False), STOPPED)

    def test_normal_exit_reads_as_done(self):
        d = self._run_dir(QUICK)
        launch(d, [sys.executable, str(d / "engine.py")])
        watch = LogWatch(d / LOG)

        def has_exit():
            watch.poll()
            return watch.exit_code is not None

        self.assertTrue(wait_for(has_exit), "no exit marker in the log")
        self.assertEqual(watch.exit_code, 0)
        self.assertEqual(watch.state(alive=False), DONE)
        self.assertIn("fake engine", (d / LOG).read_text(encoding="utf-8"))
        self.assertFalse(is_alive(d))


if __name__ == "__main__":
    unittest.main(verbosity=2)
