"""RDKit interface: geometry optimization and MMFF94 pair parameters.

This is the only module that touches RDKit, and it imports it lazily so the
rest of the GUI loads (and the test suite runs) without RDKit installed.

Two services:

- optimize_molecule(): relax a hand-built or loaded geometry with MMFF94
  (UFF fallback), keeping it planar (the engine is a 2D model) and recentring
  it on the centroid so the molecule centre sits at the origin, which is the
  frame the forcefield and the engine use (see molecule_model.h).

- mmff_pair_params(): the typed parameters needed to build a pair potential
  from the atomistic model: per-atom partial charges and the combined
  buffered-14-7 vdW parameters for every atom pair. These come straight from
  RDKit's MMFF typing. Note that RDKit's own force-field object does NOT score
  the interaction between two disjoint fragments, so the pair energy has to be
  assembled from these parameters by hand (see generate.MMFFBackend).
"""

import math
from dataclasses import dataclass

from .molecule import Atom, Molecule, bonded_pairs, connected_components


class MMFFError(Exception):
    pass


def rdkit_available() -> bool:
    """True when RDKit can be imported. The generator and the optimizer are
    disabled with a friendly note when it is missing."""
    try:
        import rdkit  # noqa: F401
        return True
    except ImportError:
        return False


def _mol_from_bonds(mol: Molecule, pairs, charge: int):
    """Build an RDKit molecule from an explicit bond list (the editor's shared
    connectivity) and let RDKit assign only the bond orders, so the perceived
    structure is identical to what the canvas draws."""
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    from rdkit.Geometry import Point3D

    rw = Chem.RWMol()
    for a in mol.atoms:
        atom = Chem.Atom(a.element)
        atom.SetNoImplicit(True)          # every atom is explicit in the model
        rw.AddAtom(atom)
    conf = Chem.Conformer(len(mol.atoms))
    for i, a in enumerate(mol.atoms):
        conf.SetAtomPosition(i, Point3D(a.x, a.y, a.z))
    for i, j in pairs:
        rw.AddBond(i, j, Chem.BondType.SINGLE)
    m = rw.GetMol()
    m.AddConformer(conf)
    rdDetermineBonds.DetermineBondOrders(m, charge=charge)
    return m


def _rdkit_mol(mol: Molecule, net_charge: int):
    """Build an RDKit molecule with perceived bond orders. Connectivity is the
    same rule the editor draws with (molecule.bonded_pairs), so the perceived
    molecule matches the canvas exactly; RDKit only fills in the bond orders."""
    if not mol.atoms:
        raise MMFFError("the molecule is empty")
    try:
        from rdkit import Chem  # noqa: F401 (import guard for a clear message)
    except ImportError:
        raise MMFFError("RDKit is not installed (pip install rdkit)")

    charge = int(net_charge)

    # reject degenerate input: real atoms are never this close (it would also be
    # drawn, and perceived, as a bond)
    for i in range(len(mol.atoms)):
        for j in range(i + 1, len(mol.atoms)):
            ai, aj = mol.atoms[i], mol.atoms[j]
            gap = math.dist((ai.x, ai.y, ai.z), (aj.x, aj.y, aj.z))
            if gap < 0.5:
                raise MMFFError(f"two atoms are almost coincident ({gap:.2f} Å "
                                "apart); check for duplicate or collapsed atoms")

    # a molecule must be one connected piece
    if len(mol.atoms) > 1 and len(connected_components(mol.atoms)) > 1:
        raise MMFFError("the molecule is not a single connected structure; "
                        "some atoms are not bonded to the rest")

    # build from the drawn bonds; RDKit only assigns the bond orders
    try:
        m = _mol_from_bonds(mol, bonded_pairs(mol.atoms), charge)
    except Exception:
        raise MMFFError("could not build a valid molecule from the bonds as "
                        "drawn; check for an atom with too many or too few "
                        "bonds, or set the right net charge for a charged species")
    if m.GetNumAtoms() != len(mol.atoms):
        raise MMFFError("atom count changed during perception")
    return m


@dataclass
class OptReport:
    method: str          # "MMFF94" or "UFF"
    converged: bool
    energy: float        # kcal/mol, final
    max_out_of_plane: float   # A, largest |z| the optimizer produced
    n_atoms: int

    def summary(self) -> str:
        state = "converged" if self.converged else "not fully converged"
        note = ""
        if self.max_out_of_plane > 0.05:
            note = (f"  ·  flattened {self.max_out_of_plane:.2f} Å out-of-plane "
                    "(the model is 2D)")
        return (f"{self.method}: {state}, energy {self.energy:.2f} kcal/mol"
                f"{note}")


def optimize_molecule(mol: Molecule, net_charge: int = 0,
                      max_iters: int = 1000) -> tuple[Molecule, OptReport]:
    """Relax the geometry with MMFF94 (UFF fallback). The result is kept planar
    (z = 0) and recentred on the centroid. Atom order and elements are
    preserved. Raises MMFFError when the molecule cannot be typed or its bonds
    cannot be perceived."""
    import numpy as np
    from rdkit.Chem import AllChem

    rdmol = _rdkit_mol(mol, net_charge)

    props = AllChem.MMFFGetMoleculeProperties(rdmol, mmffVariant="MMFF94")
    if props is not None:
        method = "MMFF94"
        ff = AllChem.MMFFGetMoleculeForceField(rdmol, props)
    elif AllChem.UFFHasAllMoleculeParams(rdmol):
        method = "UFF"
        ff = AllChem.UFFGetMoleculeForceField(rdmol)
    else:
        raise MMFFError("neither MMFF94 nor UFF can type this molecule")

    ff.Initialize()
    converged = ff.Minimize(maxIts=max_iters) == 0
    energy = float(ff.CalcEnergy())

    pos = rdmol.GetConformer().GetPositions()          # (n, 3), A
    max_z = float(np.abs(pos[:, 2]).max())
    xy = pos[:, :2] - pos[:, :2].mean(axis=0)          # centre on centroid
    atoms = [Atom(a.element, float(xy[i, 0]), float(xy[i, 1]), 0.0)
             for i, a in enumerate(mol.atoms)]
    report = OptReport(method, converged, energy, max_z, len(atoms))
    return Molecule(atoms, mol.comment), report


@dataclass
class MMFFParams:
    """Typed MMFF94 parameters for the atomistic pair potential, in the
    molecule's own frame (the origin is the rotation centre used by the
    engine). Rstar and eps are the pairwise-combined values RDKit produces,
    including the MMFF donor/acceptor contraction for hydrogen bonds."""
    xy: "object"      # (n, 2) float ndarray, A
    q: "object"       # (n,) float ndarray, e
    Rstar: "object"   # (n, n) float ndarray, A
    eps: "object"     # (n, n) float ndarray, kcal/mol
    elements: list


def mmff_pair_params(mol: Molecule, net_charge: int = 0) -> MMFFParams:
    """Extract the MMFF94 partial charges and combined vdW parameters for the
    molecule, ready for the pair-energy sum. Raises MMFFError when MMFF94 does
    not cover the molecule (there is no vdW fallback here: a potential built
    from UFF would be misleading, so we refuse instead)."""
    import numpy as np
    from rdkit.Chem import AllChem

    rdmol = _rdkit_mol(mol, net_charge)
    props = AllChem.MMFFGetMoleculeProperties(rdmol, mmffVariant="MMFF94")
    if props is None:
        raise MMFFError("MMFF94 cannot type this molecule; it likely contains "
                        "elements outside MMFF's coverage")

    n = rdmol.GetNumAtoms()
    pos = rdmol.GetConformer().GetPositions()
    xy = pos[:, :2].copy()                              # as-is: origin = centre
    q = np.array([props.GetMMFFPartialCharge(i) for i in range(n)])
    Rstar = np.empty((n, n))
    eps = np.empty((n, n))
    for i in range(n):
        for j in range(n):
            p = props.GetMMFFVdWParams(i, j)           # (R*_ij, eps_ij, ...)
            Rstar[i, j] = p[0]
            eps[i, j] = p[1]
    elements = [a.element for a in mol.atoms]
    return MMFFParams(xy=xy, q=q, Rstar=Rstar, eps=eps, elements=elements)
