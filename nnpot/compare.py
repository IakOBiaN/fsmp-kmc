import numpy as np

from fsmp_gui.forcefield import ForcefieldGrid


class ComparisonError(Exception):
    pass


def _stride(fine, coarse, what):
    ratio = coarse / fine
    steps = int(round(ratio))
    if steps < 1 or abs(ratio - steps) > 1e-6:
        raise ComparisonError(
            f"the {what} steps are not commensurate: {fine} and {coarse}")
    return steps


def _angle_points(info, period, count):
    steps = _stride(info.da, period / count, "angle")
    index = np.arange(count) * steps
    if index[-1] >= info.n_ang:
        raise ComparisonError("the angular grids do not cover the same period")
    return index


def _distance_points(info, values):
    offset = (values - info.min_dist) / info.dr
    index = np.rint(offset).astype(int)
    if np.abs(offset - index).max() > 1e-6:
        raise ComparisonError(
            f"{info.path.name} has no grid points at the shared distances")
    return index


def _slice(info, data, r_values, angle_index):
    rows = _distance_points(info, r_values)
    block = np.asarray(data[rows], dtype=float)
    return block[:, angle_index][:, :, angle_index]


def load_pair(path_a, path_b):
    a, b = ForcefieldGrid.open(path_a), ForcefieldGrid.open(path_b)
    if abs(a.info.fold - b.info.fold) > 1e-6:
        raise ComparisonError(
            f"different folding periods: {a.info.fold} and {b.info.fold}")
    period = a.info.fold

    coarse_da = max(a.info.da, b.info.da)
    count = int(round(period / coarse_da))
    if abs(count * coarse_da - period) > 1e-6:
        raise ComparisonError("the angular step does not divide the folded period")

    r_min = max(a.info.min_dist, b.info.min_dist)
    r_max = min(a.info.r_max, b.info.r_max)
    if r_max <= r_min:
        raise ComparisonError("the distance ranges do not overlap")
    dr = max(a.info.dr, b.info.dr)
    r_values = r_min + np.arange(int(np.floor((r_max - r_min) / dr + 1e-6)) + 1) * dr

    block_a = _slice(a.info, a._data, r_values, _angle_points(a.info, period, count))
    block_b = _slice(b.info, b._data, r_values, _angle_points(b.info, period, count))
    return block_a, block_b, r_values, period, count


def _turned(block, shift, mirror):
    index = np.arange(block.shape[1])
    if mirror:
        index = (-index) % block.shape[1]
    index = (index + shift) % block.shape[1]
    return block[:, index][:, :, index]


def compare(path_a, path_b, threshold=5000.0):
    block_a, block_b, r_values, period, count = load_pair(path_a, path_b)
    best = None
    for mirror in (False, True):
        for shift in range(count):
            turned = _turned(block_b, shift, mirror)
            mask = (block_a < threshold) & (turned < threshold)
            if not mask.any():
                continue
            difference = block_a[mask] - turned[mask]
            rmse = float(np.sqrt((difference ** 2).mean()))
            if best is None or rmse < best["rmse"]:
                best = {
                    "shift": shift * period / count,
                    "mirror": mirror,
                    "rmse": rmse,
                    "bias": float(difference.mean()),
                    "worst": float(np.abs(difference).max()),
                    "points": int(mask.sum()),
                }
    if best is None:
        raise ComparisonError("no points below the threshold to compare")

    turned = _turned(block_b, int(round(best["shift"] * count / period)), best["mirror"])
    for name, block in (("a", block_a), ("b", turned)):
        flat = int(np.argmin(block))
        i, j, k = np.unravel_index(flat, block.shape)
        best[name + "_min"] = float(block[i, j, k])
        best[name + "_at"] = (float(r_values[i]), j * period / count, k * period / count)
    best["r_range"] = (float(r_values[0]), float(r_values[-1]))
    best["shape"] = tuple(int(n) for n in block_a.shape)
    best["threshold"] = threshold
    return best


def report(result, name_a, name_b):
    lines = [
        f"compared on {result['shape'][0]} distances from {result['r_range'][0]:.2f} "
        f"to {result['r_range'][1]:.2f} A, {result['shape'][1]} angles each axis",
        f"best alignment: both orientations turned by {result['shift']:.1f} deg"
        + (", mirrored" if result["mirror"] else ""),
        f"points below {result['threshold'] / 1000:.0f} kJ/mol: {result['points']}",
        f"rmse  {result['rmse']:9.1f} J/mol",
        f"bias  {result['bias']:9.1f} J/mol",
        f"worst {result['worst']:9.1f} J/mol",
        f"{name_a}: minimum {result['a_min']:.1f} J/mol at r = {result['a_at'][0]:.2f} A, "
        f"a1 = {result['a_at'][1]:.1f}, a2 = {result['a_at'][2]:.1f} deg",
        f"{name_b}: minimum {result['b_min']:.1f} J/mol at r = {result['b_at'][0]:.2f} A, "
        f"a1 = {result['b_at'][1]:.1f}, a2 = {result['b_at'][2]:.1f} deg",
    ]
    return "\n".join(lines)
