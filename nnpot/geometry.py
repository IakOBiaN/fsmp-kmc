import numpy as np

ATOMIC_NUMBER = {
    "H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Si": 14, "P": 15,
    "S": 16, "Cl": 17, "As": 33, "Se": 34, "Br": 35, "I": 53,
}


def rotate(xy, angle_deg):
    a = np.deg2rad(angle_deg)
    c, s = np.cos(a), np.sin(a)
    x, y = xy[:, 0], xy[:, 1]
    return np.stack([x * c - y * s, x * s + y * c], axis=1)


def rotate_batch(xy, angles_rad):
    c, s = np.cos(angles_rad), np.sin(angles_rad)
    x, y = xy[:, 0], xy[:, 1]
    return np.stack([np.outer(c, x) - np.outer(s, y),
                     np.outer(s, x) + np.outer(c, y)], axis=-1)


def _match(xy, elements, angle_deg):
    rotated = rotate(xy, angle_deg)
    n = len(xy)
    perm = np.full(n, -1, dtype=int)
    worst = 0.0
    taken = set()
    for i in range(n):
        best, best_d = -1, np.inf
        for j in range(n):
            if j in taken or elements[j] != elements[i]:
                continue
            d = float(np.hypot(*(rotated[i] - xy[j])))
            if d < best_d:
                best, best_d = j, d
        if best < 0:
            return None, np.inf
        perm[i] = best
        taken.add(best)
        worst = max(worst, best_d)
    return perm, worst


def symmetry_error(xy, elements, order):
    if order <= 1:
        return 0.0
    worst = 0.0
    for k in range(1, order):
        _, d = _match(np.asarray(xy, float), list(elements), 360.0 * k / order)
        worst = max(worst, d)
    return worst


def rotational_order(xy, elements, tol=0.05, max_order=12):
    xy = np.asarray(xy, float)
    elements = list(elements)
    for n in range(max_order, 1, -1):
        if symmetry_error(xy, elements, n) <= tol:
            return n
    return 1


def symmetrize(xy, elements, order):
    xy = np.asarray(xy, float)
    if order <= 1:
        return xy.copy()
    perm, _ = _match(xy, list(elements), 360.0 / order)
    if perm is None:
        raise ValueError(f"the molecule cannot be matched onto itself under C{order}")
    total = xy.copy()
    current = np.arange(len(xy))
    for k in range(1, order):
        current = perm[current]
        total = total + rotate(xy[current], -360.0 * k / order)
    return total / order


def orbit_map(p, da):
    i, j = np.meshgrid(np.arange(p), np.arange(p), indexing="ij")
    flat = i * p + j
    half = 180.0 / da
    if abs(half - round(half)) > 1e-9:
        canon = flat
    else:
        shift = int(round(half)) % p
        partner = ((j + shift) % p) * p + ((i + shift) % p)
        canon = np.minimum(flat, partner)
    uniq, inverse = np.unique(canon.ravel(), return_inverse=True)
    return uniq, inverse.reshape(p, p)


def pair_distances(xy, r, t1, t2):
    a = rotate_batch(xy, t1)
    b = rotate_batch(xy, t2) + np.array([r, 0.0])
    d = a[:, :, None, :] - b[:, None, :, :]
    return np.sqrt(np.einsum("kijc,kijc->kij", d, d))


def dimer_coords(xy, r, t1, t2):
    a = rotate_batch(xy, t1)
    b = rotate_batch(xy, t2) + np.array([r, 0.0])
    flat = np.concatenate([a, b], axis=1)
    return np.concatenate([flat, np.zeros(flat.shape[:2] + (1,))], axis=2)
