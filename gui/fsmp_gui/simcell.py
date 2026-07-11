"""The elongated simulation cell and its fields, mirrored from the engine
(generate_elongated_cell in StructureGenerator.h and fields.h) so the
Simulation cell tab can preview exactly what a run would build: a crystal
slab of uc_in_x x uc_in_y unit cells centered in a periodic box with free
space on both sides, the damping field lambda(x) and the external potential
u_ext(x).

Both fields depend only on the reduced distance from the cell center,
ksi = 32|x - Lx/2|/Lx (16 at the cell edge): full interactions below 5, a
smooth ramp to the transition-zone plateau at 7..9 (lambda = T/T_tz,
u_ext = u_m), a second ramp to the gas zone above 11 (lambda = lambdam^2,
u_ext = u_m). Where the crystal slab edge lands among these zones is set by
free_space, which is the whole point of previewing them together.
"""

import math
from dataclasses import dataclass

# zone boundaries in ksi units (fields.h); the cell edge is at 16
KSI_CRYSTAL = 5.0
KSI_PLATEAU = 7.0
KSI_PLATEAU_END = 9.0
KSI_GAS = 11.0
KSI_EDGE = 16.0

ENGINE_BUFFER = 5000   # vector<state> coordinates in program_body.cpp


@dataclass
class SimCellSpec:
    """The parameters that shape the elongated cell (engine key names)."""
    uc_in_x: int = 22
    uc_in_y: int = 6
    free_space: float = 0.24        # fraction of Lx free on each side
    lambdam: float = 0.0            # sqrt(lambda) in the gas zone, as in the engine
    temp_transition: float = 900.0  # temperature_in_transition_zone, K


def ksi_at(x: float, lx: float) -> float:
    return 32.0 * abs(x - lx / 2.0) / lx


def damping_sqrt(ksi: float, lambda0: float, lambdam: float) -> float:
    """sqrt(lambda)(ksi): the exact piecewise form of fields.h."""
    if ksi < KSI_CRYSTAL:
        return 1.0
    if ksi < KSI_PLATEAU:
        return lambda0 + (1.0 - lambda0) * (7.0 - ksi) ** 2 * (ksi - 4.0) / 4.0
    if ksi < KSI_PLATEAU_END:
        return lambda0
    if ksi < KSI_GAS:
        return lambdam + (lambda0 - lambdam) * (11.0 - ksi) ** 2 * (ksi - 8.0) / 4.0
    return lambdam


def external_u(ksi: float, um: float) -> float:
    """u_ext(ksi) in the units of um: the exact piecewise form of fields.h."""
    if ksi < KSI_CRYSTAL:
        return 0.0
    if ksi > KSI_PLATEAU:
        return um
    return um * (1.0 - (7.0 - ksi) ** 2 * (ksi - 4.0) / 4.0)


def lambda_profile(x: float, lx: float, temp: float, spec: SimCellSpec) -> float:
    """The physical damping lambda(x): 1 in the crystal, T/T_tz on the
    plateau, lambdam^2 in the gas. The engine stores square roots
    (lambda0 = sqrt(T/T_tz)); the physical coefficient is their square."""
    lambda0 = math.sqrt(temp / spec.temp_transition)
    return damping_sqrt(ksi_at(x, lx), lambda0, spec.lambdam) ** 2


def u_profile(x: float, lx: float, um: float) -> float:
    """u_ext(x) in the units of um (J/mol in the parameter file)."""
    return external_u(ksi_at(x, lx), um)


class SimCellLayout:
    """Geometry of the elongated cell built from a cell_x x cell_y unit cell
    holding n_in_cell molecules (generate_elongated_cell)."""

    def __init__(self, cell_x: float, cell_y: float, n_in_cell: int,
                 spec: SimCellSpec):
        self.spec = spec
        self.slab_x = cell_x * spec.uc_in_x
        self.lx = self.slab_x / (1.0 - 2.0 * spec.free_space)
        self.ly = cell_y * spec.uc_in_y
        self.slab_x0 = (self.lx - self.slab_x) / 2.0
        self.n_molecules = n_in_cell * spec.uc_in_x * spec.uc_in_y

    def zone_x(self, ksi: float) -> tuple[float, float]:
        """The two x positions at a given ksi, symmetric about the center."""
        d = ksi * self.lx / 32.0
        return self.lx / 2.0 - d, self.lx / 2.0 + d


def tile_unit_cell(unit_cell: dict, spec: SimCellSpec):
    """(layout, [(x, y, phi), ...]): the saved unit cell tiled
    uc_in_x x uc_in_y times and centered in the elongated box. The engine
    centers the slab by its center of mass; centering the slab box is the
    same up to a sub-cell offset and is enough for the preview."""
    mols = unit_cell["molecules"]
    cx, cy = unit_cell["cell_x"], unit_cell["cell_y"]
    layout = SimCellLayout(cx, cy, len(mols), spec)
    placements = []
    for i in range(spec.uc_in_x):
        for j in range(spec.uc_in_y):
            for m in mols:
                placements.append((layout.slab_x0 + i * cx + m["x"] % cx,
                                   j * cy + m["y"] % cy, m["phi"]))
    return layout, placements
