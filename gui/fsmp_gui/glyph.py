"""Turn the project molecule model into a drawable glyph: a list of coloured
disks (x, y, color, radius) in the molecule's local frame. Shared by the
potential viewer and the unit-cell editor."""

from . import theme
from .elements import covalent_radius, element_color
from .molecule import Molecule
from .project import Project


def model_glyph(project: Project):
    """Glyph of the project's atomistic model, or None when the model is
    absent or unreadable (callers draw fallback_glyph then)."""
    entry = project.atomistic
    if entry is not None:
        try:
            mol = Molecule.load_xyz(project.model_path(entry))
            return [(a.x, a.y, element_color(a.element),
                     min(max(0.4 * covalent_radius(a.element) + 0.15, 0.3), 0.8))
                    for a in mol.atoms]
        except (OSError, ValueError):
            pass
    return None


def fallback_glyph():
    return [(0.0, 0.0, theme.ACCENT, 0.5)]
