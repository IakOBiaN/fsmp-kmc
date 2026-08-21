import numpy as np

EV_TO_JMOL = 96485.33212


class AIMNet2:

    name = "AIMNet2"

    def __init__(self, model="aimnet2", coulomb="simple", charge=0.0):
        from aimnet2calc import AIMNet2Calculator
        self.calc = AIMNet2Calculator(model)
        self.model_name = model
        self.charge = float(charge)
        self.device = str(self.calc.device)
        self.short_cutoff = float(self.calc.cutoff)
        if not self.calc.lr:
            self.coulomb = None
            self.cutoff = self.short_cutoff
        else:
            self.calc.set_lrcoulomb_method(coulomb)
            self.coulomb = coulomb
            self.cutoff = None if coulomb == "simple" else float(self.calc.cutoff_lr)

    def describe(self):
        return {
            "model": self.model_name,
            "device": self.device,
            "coulomb": self.coulomb,
            "short_cutoff": self.short_cutoff,
            "declared_cutoff": self.cutoff,
        }

    def energies(self, z, coords):
        coords = np.ascontiguousarray(coords, dtype=np.float32)
        if coords.ndim == 2:
            coords = coords[None]
        k, n, _ = coords.shape
        numbers = np.broadcast_to(np.asarray(z, dtype=np.int64), (k, n))
        out = self.calc.eval({
            "coord": coords,
            "numbers": np.ascontiguousarray(numbers),
            "charge": np.full(k, self.charge, dtype=np.float32),
        })
        energy = out["energy"].detach().cpu().numpy().reshape(-1).astype(float)
        return energy * EV_TO_JMOL
