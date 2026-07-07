"""RDKit-backed chemistry: classical force-field geometry optimization
(MMFF94 / UFF) and partial charges (Gasteiger / MMFF94).

RDKit is a listed dependency, but everything degrades gracefully when it
is missing so that the rest of the GUI stays usable.
"""

from dataclasses import dataclass

from .molecule import Atom, Molecule

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDetermineBonds, rdPartialCharges
    from rdkit.Geometry import Point3D
    HAVE_RDKIT = True
except ImportError:
    HAVE_RDKIT = False

KCAL_TO_KJ = 4.184

OPTIMIZE_METHODS = ["MMFF94", "UFF"]
CHARGE_METHODS = ["Gasteiger", "MMFF94"]


class ChemError(Exception):
    pass


@dataclass
class OptimizeResult:
    molecule: Molecule
    e_before: float  # kJ/mol
    e_after: float   # kJ/mol
    converged: bool


def _to_rdkit(molecule: Molecule, charge: int):
    if not HAVE_RDKIT:
        raise ChemError("RDKit is not installed (pip install rdkit)")
    if not molecule.atoms:
        raise ChemError("the molecule is empty")
    rw = Chem.RWMol()
    conf = Chem.Conformer(len(molecule.atoms))
    for i, a in enumerate(molecule.atoms):
        try:
            idx = rw.AddAtom(Chem.Atom(a.element))
        except Exception:
            raise ChemError(f"unknown element '{a.element}' (atom {i + 1})")
        conf.SetAtomPosition(idx, Point3D(a.x, a.y, a.z))
    mol = rw.GetMol()
    mol.AddConformer(conf)
    try:
        rdDetermineBonds.DetermineBonds(mol, charge=charge)
    except (ValueError, RuntimeError) as e:
        raise ChemError("could not determine bonds from the geometry: "
                        f"{e}\nCheck bond lengths, missing hydrogens and "
                        "the total charge.")
    return mol


def _from_rdkit(mol, template: Molecule) -> Molecule:
    conf = mol.GetConformer()
    atoms = []
    for i, a in enumerate(template.atoms):
        p = conf.GetAtomPosition(i)
        atoms.append(Atom(a.element, round(p.x, 6), round(p.y, 6), round(p.z, 6)))
    return Molecule(atoms, template.comment)


def _force_field(mol, method: str):
    if method == "MMFF94":
        props = AllChem.MMFFGetMoleculeProperties(mol)
        if props is None:
            raise ChemError("MMFF94 has no parameters for this molecule; try UFF")
        ff = AllChem.MMFFGetMoleculeForceField(mol, props)
    elif method == "UFF":
        ff = AllChem.UFFGetMoleculeForceField(mol)
    else:
        raise ChemError(f"unknown method '{method}'")
    if ff is None:
        raise ChemError(f"{method} setup failed for this molecule")
    return ff


def optimize(molecule: Molecule, method: str = "MMFF94", charge: int = 0,
             max_steps: int = 500) -> OptimizeResult:
    mol = _to_rdkit(molecule, charge)
    ff = _force_field(mol, method)
    e_before = ff.CalcEnergy() * KCAL_TO_KJ
    code = ff.Minimize(maxIts=max_steps)
    e_after = ff.CalcEnergy() * KCAL_TO_KJ
    return OptimizeResult(_from_rdkit(mol, molecule), e_before, e_after,
                          code == 0)


def partial_charges(molecule: Molecule, method: str = "Gasteiger",
                    charge: int = 0) -> list[float]:
    mol = _to_rdkit(molecule, charge)
    if method == "Gasteiger":
        rdPartialCharges.ComputeGasteigerCharges(mol)
        values = [a.GetDoubleProp("_GasteigerCharge") for a in mol.GetAtoms()]
    elif method == "MMFF94":
        props = AllChem.MMFFGetMoleculeProperties(mol)
        if props is None:
            raise ChemError("MMFF94 has no parameters for this molecule")
        values = [props.GetMMFFPartialCharge(i) for i in range(mol.GetNumAtoms())]
    else:
        raise ChemError(f"unknown method '{method}'")
    if any(v != v for v in values):  # NaN guard
        raise ChemError("the charge calculation produced NaN values")
    return [round(v, 5) for v in values]
