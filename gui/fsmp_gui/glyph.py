"""Turn the project molecule model into a drawable glyph: a list of coloured
disks (x, y, color, radius) in the molecule's local frame. Shared by the
potential viewer and the unit-cell editor."""

from . import theme
from .canvas import site_color
from .elements import covalent_radius, element_color
from .molecule import Molecule
from .project import Project
from .sitemodel import SiteModel


def model_glyph(project: Project):
    """Glyph for visualization: prefer the atomistic model (it looks like the
    real molecule); fall back to the site model. None if neither is readable."""
    entry = project.atomistic
    if entry is not None:
        try:
            mol = Molecule.load_xyz(project.model_path(entry))
            return [(a.x, a.y, element_color(a.element),
                     min(max(0.4 * covalent_radius(a.element) + 0.15, 0.3), 0.8))
                    for a in mol.atoms]
        except (OSError, ValueError):
            pass
    entry = project.site
    if entry is not None:
        try:
            sm = SiteModel.load(project.model_path(entry))
            return [(s.x, s.y, site_color(s.q).name(),
                     0.32 if s.sigma <= 0 else min(max(0.2 * s.sigma + 0.1, 0.32), 0.9))
                    for s in sm.sites]
        except (OSError, ValueError):
            pass
    return None


def fallback_glyph():
    return [(0.0, 0.0, theme.ACCENT, 0.5)]
