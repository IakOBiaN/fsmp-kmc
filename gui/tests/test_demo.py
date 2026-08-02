"""The one-click demonstration: installing it must produce a self-contained,
runnable project copy, and the start page must lead straight to the Run tab.

Run from the repository root:
    gui/.venv/Scripts/python gui/tests/test_demo.py
"""

import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from PySide6.QtWidgets import QApplication

from fsmp_gui import demo
from fsmp_gui.forcefield import read_header
from fsmp_gui.project import Project

_app = QApplication.instance() or QApplication([])


@unittest.skipUnless(demo.available(), "demonstration data not present")
class TestInstall(unittest.TestCase):
    def test_installed_copy_is_self_contained(self):
        with tempfile.TemporaryDirectory() as td:
            project, created = demo.install(Path(td))
            self.assertTrue(created)
            self.assertEqual(project.root, Path(td) / demo.PROJECT_NAME)

            # the potential travels with the copy, stored relative to it
            self.assertFalse(Path(project.potential["path"]).is_absolute())
            grid = project.potential_path()
            self.assertTrue(grid.is_file())
            self.assertEqual(grid.parent, project.root / demo.POTENTIAL_DIR)
            self.assertEqual(read_header(grid).version, 2)

            # and everything a run needs came along
            self.assertIsNotNone(project.atomistic)
            self.assertTrue(project.model_path(project.atomistic).is_file())
            self.assertTrue(project.unit_cell["molecules"])
            self.assertIsNotNone(project.simulation_cell)
            self.assertIsNotNone(project.simulation)

    def test_moving_the_copy_keeps_it_working(self):
        """No path in the installed manifest may point outside the folder."""
        with tempfile.TemporaryDirectory() as td:
            demo.install(Path(td) / "here")
            moved = Path(td) / "elsewhere" / "renamed"
            moved.parent.mkdir()
            shutil.move(str(Path(td) / "here" / demo.PROJECT_NAME), str(moved))
            project = Project.open(moved)
            self.assertTrue(project.potential_path().is_file())

    def test_second_install_reopens_instead_of_duplicating(self):
        with tempfile.TemporaryDirectory() as td:
            first, created_first = demo.install(Path(td))
            (first.root / "runs").mkdir()          # a run the user started
            second, created_second = demo.install(Path(td))
            self.assertTrue(created_first)
            self.assertFalse(created_second)
            self.assertEqual(first.root, second.root)
            self.assertTrue((second.root / "runs").is_dir())
            self.assertEqual(len(list(Path(td).iterdir())), 1)

    def test_refuses_a_foreign_folder(self):
        with tempfile.TemporaryDirectory() as td:
            occupied = Path(td) / demo.PROJECT_NAME
            occupied.mkdir()
            (occupied / "notes.txt").write_text("mine")
            with self.assertRaises(demo.DemoError):
                demo.install(Path(td))


@unittest.skipUnless(demo.available(), "demonstration data not present")
class TestStartPage(unittest.TestCase):
    def test_button_opens_the_demo_on_the_run_tab(self):
        from fsmp_gui.main_window import MainWindow

        # ignore_cleanup_errors: the Potentials tab memory-maps the grid to
        # draw it, and Windows keeps the mapping until the widget is gone
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
            window = MainWindow()
            window.demo_location = lambda: Path(td)
            self.assertTrue(window.start_page.demo_btn.isEnabled())
            window.start_page.demo_btn.click()

            view = window.project_view
            self.assertIsNotNone(view, "the demonstration did not open")
            self.assertIs(view.tabs.currentWidget(), view.run_tab)
            self.assertEqual(view.project.root, Path(td) / demo.PROJECT_NAME)
            # the Run tab is ready: nothing left to fetch or fill in
            view.run_tab._check_prereqs()
            self.assertEqual(view.run_tab.prereq.text(), "")

            window.close_project()
            window.close()
            _app.processEvents()


if __name__ == "__main__":
    unittest.main(verbosity=2)
