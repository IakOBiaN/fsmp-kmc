"""The generation worker lifecycle: every run ends in exactly one of
finished_ok / failed / cancelled, and the page's controls come back after a
cancel. This locks the cancel path, which used to emit nothing and leave the
page stuck on "Cancelling…" with Generate disabled until a restart.

Run from the repository root:
    gui/.venv/Scripts/python gui/tests/test_create_potential.py
"""

import os
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from fsmp_gui.forcefield import read_header
from fsmp_gui.generate import GridSpec, generate
from fsmp_gui.project import Project
from fsmp_gui.tabs.create_potential_tab import AtomisticPage, GenerateWorker

_app = QApplication.instance() or QApplication([])

SPEC = GridSpec(5.0, 7.0, 1.0, 30.0, 360.0, use_float=True)   # 3 tiny rows


class GateBackend:
    """First slab call signals `entered` and blocks until `release`, so a test
    can cancel at a known point instead of racing the worker."""

    def __init__(self):
        self.entered = threading.Event()
        self.release = threading.Event()

    def slab(self, r, ang):
        self.entered.set()
        self.release.wait(10.0)
        return np.zeros((len(ang), len(ang)))


class FreeBackend:
    def slab(self, r, ang):
        return np.zeros((len(ang), len(ang)))


def collect(worker):
    """Connect all outcome signals to plain lists. Direct connections run in
    the worker thread, so no event loop is needed to observe them."""
    got = {"ok": [], "failed": [], "cancelled": []}
    worker.finished_ok.connect(got["ok"].append, Qt.DirectConnection)
    worker.failed.connect(got["failed"].append, Qt.DirectConnection)
    worker.cancelled.connect(lambda: got["cancelled"].append(True),
                             Qt.DirectConnection)
    return got


class TestGenerateWorker(unittest.TestCase):
    def test_cancel_reports_cancelled_and_removes_the_file(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
            out = Path(td) / "toy.v2.bin"
            gate = GateBackend()
            worker = GenerateWorker(lambda: gate, SPEC, out)
            got = collect(worker)
            worker.start()
            self.assertTrue(gate.entered.wait(10.0), "backend never entered")
            worker.cancel()
            gate.release.set()
            self.assertTrue(worker.wait(10000), "worker did not finish")
            self.assertEqual(got, {"ok": [], "failed": [], "cancelled": [True]})
            self.assertFalse(out.exists(), "partial file left behind")

    def test_completion_reports_finished_ok(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
            out = Path(td) / "toy.v2.bin"
            worker = GenerateWorker(FreeBackend, SPEC, out)
            got = collect(worker)
            worker.start()
            self.assertTrue(worker.wait(10000), "worker did not finish")
            self.assertEqual(got, {"ok": [str(out)], "failed": [],
                                   "cancelled": []})
            self.assertEqual(read_header(out).n_dist, SPEC.n_dist)

    def test_backend_error_reports_failed(self):
        def broken():
            raise ValueError("no parameters for this molecule")

        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
            worker = GenerateWorker(broken, SPEC, Path(td) / "toy.v2.bin")
            got = collect(worker)
            worker.start()
            self.assertTrue(worker.wait(10000), "worker did not finish")
            self.assertEqual(got["failed"], ["no parameters for this molecule"])
            self.assertEqual((got["ok"], got["cancelled"]), ([], []))

    def test_generate_reports_completion(self):
        """The generate() contract the worker builds on: True for a written
        grid, False (and no file) when cancelled."""
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
            out = Path(td) / "toy.v2.bin"
            self.assertTrue(generate(FreeBackend(), SPEC, out))
            self.assertTrue(out.exists())
            self.assertFalse(generate(FreeBackend(), SPEC, out,
                                      cancel=lambda: True))
            self.assertFalse(out.exists())

    def test_cancel_after_the_last_row_is_a_completed_grid(self):
        """A cancel that lands after generate() wrote everything must read as
        success: the file is whole, pretending otherwise would discard it."""
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
            out = Path(td) / "toy.v2.bin"
            holder = {}

            class LastRowCancel(FreeBackend):
                def __init__(self):
                    self.calls = 0

                def slab(self, r, ang):
                    self.calls += 1
                    if self.calls == SPEC.n_dist:   # during the final row
                        holder["worker"].cancel()
                    return super().slab(r, ang)

            worker = GenerateWorker(LastRowCancel, SPEC, out)
            holder["worker"] = worker
            got = collect(worker)
            worker.run()             # synchronously: the outcome decision only
            self.assertEqual(got["ok"], [str(out)])
            self.assertEqual(got["cancelled"], [])
            self.assertTrue(out.exists())


class TestPageWiring(unittest.TestCase):
    """Drive the real AtomisticPage through generate -> cancel and require the
    controls to come back, exercising the actual signal wiring."""

    def test_cancel_restores_the_page(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
            project = Project.create(Path(td) / "proj", "t")
            page = AtomisticPage(project)
            out = Path(td) / "toy.v2.bin"
            gate = GateBackend()
            page._prepare_generation = lambda: ((lambda: gate), out)

            page._generate()
            worker = page._worker
            self.assertIsNotNone(worker)
            self.assertFalse(page.gen_btn.isEnabled())
            self.assertTrue(page.cancel_btn.isEnabled())
            self.assertTrue(gate.entered.wait(10.0), "backend never entered")

            page._cancel()
            self.assertFalse(page.cancel_btn.isEnabled())
            gate.release.set()

            deadline = time.monotonic() + 10.0
            while page._worker is not None and time.monotonic() < deadline:
                _app.processEvents()
                time.sleep(0.01)

            self.assertIsNone(page._worker, "page never learned of the cancel")
            self.assertTrue(page.gen_btn.isEnabled())
            self.assertFalse(page.cancel_btn.isEnabled())
            self.assertEqual(page.gen_status.text(), "Cancelled.")
            self.assertFalse(out.exists())
            self.assertTrue(worker.wait(10000))


if __name__ == "__main__":
    unittest.main(verbosity=2)
