"""The Run tab's 'Save settings' button persists the form into the project,
so a prepared project reopens with its intended parameters (the basis for a
folder of ready-to-run example projects). Runs under the offscreen platform.

Run from the repository root:
    gui/.venv/Scripts/python gui/tests/test_run_tab.py
"""

import os
import sys
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from PySide6.QtWidgets import QApplication

from fsmp_gui.project import Project
from fsmp_gui.tabs.run_tab import RunTab

_app = QApplication.instance() or QApplication([])


class TestSaveSettings(unittest.TestCase):
    def test_save_button_persists_form_to_disk(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td) / "proj"
            project = Project.create(root, "p")
            tab = RunTab(project)
            self.assertIsNone(project.simulation)   # nothing saved yet

            tab.temp_from.setValue(250.0)
            tab.temp_to.setValue(500.0)
            tab.nsteps.setValue(1234567)
            tab.kmc.setChecked(False)
            tab.sigma_mode.setCurrentText("min_dist")
            tab.mask_enable.setChecked(True)
            tab.seed.setValue(42)
            tab._save_settings()

            form = Project.open(root).simulation   # reload from project.json
            self.assertIsNotNone(form)
            self.assertEqual(form["temp_from"], 250.0)
            self.assertEqual(form["temp_to"], 500.0)
            self.assertEqual(form["nSteps"], 1234567)
            self.assertFalse(form["kMC"])
            self.assertEqual(form["sigma_mode"], "min_dist")
            self.assertTrue(form["mask"])
            self.assertEqual(form["seed"], 42)

    def test_reopened_project_loads_the_saved_form(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td) / "proj"
            project = Project.create(root, "p")
            tab = RunTab(project)
            tab.temp_from.setValue(275.0)
            tab.delta.setValue(3.5)
            tab.widom.setChecked(True)
            tab.mask_penalty.setValue(30000.0)
            tab._save_settings()

            # a fresh tab on the reopened project restores the widgets, which
            # is what the operator sees when opening a prepared project
            tab2 = RunTab(Project.open(root))
            tab2._load_form_defaults()
            self.assertEqual(tab2.temp_from.value(), 275.0)
            self.assertEqual(tab2.delta.value(), 3.5)
            self.assertTrue(tab2.widom.isChecked())
            self.assertEqual(tab2.mask_penalty.value(), 30000.0)

    def test_editing_clears_the_saved_confirmation(self):
        with tempfile.TemporaryDirectory() as td:
            project = Project.create(Path(td) / "proj", "p")
            tab = RunTab(project)
            tab._save_settings()
            self.assertIn("saved", tab.save_status.text().lower())
            tab.temp_from.setValue(tab.temp_from.value() + 10)
            self.assertEqual(tab.save_status.text().strip(), "")


class TestPrereqs(unittest.TestCase):
    """The Run tab's readiness banner. A sample project ships with a potential
    attached but the heavy file is a separate download, so 'attached' must not
    be mistaken for 'the file is here'."""

    def test_unattached_potential_is_flagged(self):
        with tempfile.TemporaryDirectory() as td:
            project = Project.create(Path(td) / "proj", "p")
            tab = RunTab(project)
            tab._check_prereqs()
            self.assertIn("potential (tab 3)", tab.prereq.text())

    def test_attached_but_missing_file_is_flagged(self):
        with tempfile.TemporaryDirectory() as td:
            project = Project.create(Path(td) / "proj", "p")
            # attach a potential whose file is not present (the download case)
            project.set_potential("ghost", Path(td) / "ghost.v2.bin")
            tab = RunTab(project)
            tab._check_prereqs()
            text = tab.prereq.text()
            self.assertIn("forcefields", text)          # tells the user what to do
            self.assertNotIn("potential (tab 3)", text)  # it IS attached

    def test_present_file_is_not_flagged(self):
        with tempfile.TemporaryDirectory() as td:
            project = Project.create(Path(td) / "proj", "p")
            grid = Path(td) / "here.v2.bin"
            grid.write_bytes(b"\x00" * 64)
            project.set_potential("here", grid)
            tab = RunTab(project)
            tab._check_prereqs()
            self.assertNotIn("potential", tab.prereq.text())


if __name__ == "__main__":
    unittest.main(verbosity=2)
