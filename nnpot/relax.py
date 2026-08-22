import numpy as np


def relax_planar(elements, xy, calculator, fmax=0.01, steps=500, logfile=None):
    from ase import Atoms
    from ase.constraints import FixedPlane
    from ase.optimize import BFGS

    elements = list(elements)
    positions = np.concatenate([np.asarray(xy, dtype=float),
                                np.zeros((len(elements), 1))], axis=1)
    atoms = Atoms(symbols=elements, positions=positions)
    atoms.set_constraint(FixedPlane(list(range(len(elements))), (0.0, 0.0, 1.0)))
    atoms.calc = calculator.ase()
    optimizer = BFGS(atoms, logfile=logfile)
    converged = bool(optimizer.run(fmax=fmax, steps=steps))
    out = atoms.get_positions()[:, :2]
    out = out - out.mean(axis=0)
    return out, converged, int(optimizer.get_number_of_steps())
