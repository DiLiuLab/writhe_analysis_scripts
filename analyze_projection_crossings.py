#!/usr/bin/env python3
"""
Analyze projection-dependent crossings and writhe of a closed 3D curve.

Input:
    XYZ-like text file containing x y z coordinates of a closed curve.
Output:
    A CSV table with one row per viewing direction and a text summary.
    Output names are derived from the input basename by adding a suffix.

The script distinguishes:
    directional signed crossing sum w(n) = sum of crossing signs,
    unsigned crossing count N(n), and
    their averages over many directions.  The spherical average of w(n)
    approaches the 3D writhe Wr, whereas the average of N(n) is the average
    crossing number (ACN), a different quantity.

Example:
    python analyze_projection_crossings.py examples/sc_W3_ALocal.xyz --n-directions 10000 \
        --analysis-points 1000 --include-axis-views

Run without arguments, or with --gui, to enter the parameters in a window.

Python 3.9+; requires NumPy.  Shapely 2.x is recommended for speed.
"""

import argparse
import csv
import math
import os
import sys
from collections import Counter
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np

from writhe_gui import add_gui_argument, gui_requested, launch_gui

try:
    import shapely
    from shapely import STRtree
    HAVE_SHAPELY = True
except Exception:
    HAVE_SHAPELY = False


def add_suffix(path: str, suffix: str, ext: Optional[str] = None) -> str:
    """Add suffix before the extension of path."""
    root, old_ext = os.path.splitext(path)
    return root + suffix + (old_ext if ext is None else ext)


def read_xyz(path: str) -> np.ndarray:
    """Read the first three numeric fields from each non-comment line."""
    pts: List[List[float]] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("!"):
                continue
            vals: List[float] = []
            for token in line.split():
                try:
                    vals.append(float(token))
                except ValueError:
                    continue
            if len(vals) >= 3:
                pts.append(vals[:3])
    if len(pts) < 4:
        raise ValueError("Need at least four 3D points for a closed curve.")
    arr = np.asarray(pts, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise ValueError("Input contains NaN or infinite coordinates.")
    return arr


def strip_closing_duplicate(points: np.ndarray, tol: float = 1.0e-9) -> np.ndarray:
    """Remove a duplicated final point if it matches the first point."""
    pts = np.asarray(points, dtype=float)
    if len(pts) > 1 and np.linalg.norm(pts[-1] - pts[0]) <= tol:
        return pts[:-1].copy()
    return pts.copy()


def resample_closed_polyline(points: np.ndarray, n_points: int) -> np.ndarray:
    """Resample a closed polyline to n_points equally spaced by arc length."""
    pts = strip_closing_duplicate(points)
    if n_points <= 0 or n_points == len(pts):
        return pts.copy()
    if n_points < 4:
        raise ValueError("analysis-points must be at least 4.")

    nxt = np.roll(pts, -1, axis=0)
    seg = nxt - pts
    lengths = np.linalg.norm(seg, axis=1)
    total = float(lengths.sum())
    if total <= 1.0e-12:
        raise ValueError("Curve has negligible total length.")

    cum = np.concatenate(([0.0], np.cumsum(lengths)))
    q = np.linspace(0.0, total, n_points, endpoint=False)
    idx = np.searchsorted(cum, q, side="right") - 1
    idx = np.clip(idx, 0, len(pts) - 1)
    local = (q - cum[idx]) / np.maximum(lengths[idx], 1.0e-15)
    return pts[idx] + local[:, None] * seg[idx]


def fibonacci_sphere(n: int) -> np.ndarray:
    """Return approximately uniform unit vectors on the full sphere."""
    if n < 1:
        raise ValueError("n-directions must be positive.")
    k = np.arange(n, dtype=float)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    z = 1.0 - 2.0 * (k + 0.5) / n
    r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    phi = golden * k
    return np.column_stack((r * np.cos(phi), r * np.sin(phi), z))


def normalize_vector(v: Sequence[float]) -> np.ndarray:
    """Normalize a 3-vector."""
    arr = np.asarray(v, dtype=float)
    norm = float(np.linalg.norm(arr))
    if norm <= 1.0e-15:
        raise ValueError("Viewing direction cannot be the zero vector.")
    return arr / norm


def projection_basis(normal: Sequence[float]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return right-handed orthonormal screen basis (e1, e2, n)."""
    n = normalize_vector(normal)
    ref = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(ref, n))) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    e1 = ref - np.dot(ref, n) * n
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    e2 /= np.linalg.norm(e2)
    return e1, e2, n


def project_points(points: np.ndarray, normal: Sequence[float]) -> np.ndarray:
    """Orthographically project 3D points into a 2D screen basis."""
    e1, e2, _n = projection_basis(normal)
    return np.column_stack((points @ e1, points @ e2))


def cross2(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """2D scalar cross product for arrays ending in dimension 2."""
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def adjacent_segment_pairs(i: np.ndarray, j: np.ndarray, nseg: int) -> np.ndarray:
    """Boolean mask identifying identical or adjacent segments of a closed curve."""
    diff = np.abs(i - j)
    return (diff == 0) | (diff == 1) | (diff == nseg - 1)


def candidate_pairs_shapely(projected: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Use Shapely STRtree to find pairs whose open interiors cross."""
    n = len(projected)
    q = np.roll(projected, -1, axis=0)
    coords = np.stack((projected, q), axis=1)
    lines = shapely.linestrings(coords)
    tree = STRtree(lines)
    pairs = tree.query(lines, predicate="crosses")
    if pairs.size == 0:
        return np.empty(0, dtype=int), np.empty(0, dtype=int)
    i = np.asarray(pairs[0], dtype=int)
    j = np.asarray(pairs[1], dtype=int)
    keep = i < j
    i = i[keep]
    j = j[keep]
    keep = ~adjacent_segment_pairs(i, j, n)
    return i[keep], j[keep]


def candidate_pairs_numpy(projected: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Portable O(N^2) fallback for candidate pairs using bounding boxes."""
    n = len(projected)
    i, j = np.triu_indices(n, 1)
    keep = ~adjacent_segment_pairs(i, j, n)
    i = i[keep]
    j = j[keep]
    q = np.roll(projected, -1, axis=0)
    pmin = np.minimum(projected, q)
    pmax = np.maximum(projected, q)
    bbox = (
        (pmax[i, 0] >= pmin[j, 0])
        & (pmax[j, 0] >= pmin[i, 0])
        & (pmax[i, 1] >= pmin[j, 1])
        & (pmax[j, 1] >= pmin[i, 1])
    )
    return i[bbox], j[bbox]


def projection_crossing_details(
    points: np.ndarray,
    normal: Sequence[float],
    engine: str = "auto",
    eps: float = 1.0e-10,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Return detailed regular crossings for one projection.

    Returns:
        crossing_xy : (K,2) projected intersection coordinates
        signs       : (K,) crossing signs using the report convention
        seg_i       : (K,) first segment indices
        seg_j       : (K,) second segment indices
        depth_delta : (K,) signed separation (r_i-r_j).n at the crossing
    """
    pts = strip_closing_duplicate(points)
    proj = project_points(pts, normal)
    nvec = normalize_vector(normal)

    use_shapely = engine == "shapely" or (engine == "auto" and HAVE_SHAPELY)
    if use_shapely:
        i, j = candidate_pairs_shapely(proj)
    elif engine in ("auto", "numpy"):
        i, j = candidate_pairs_numpy(proj)
    else:
        raise ValueError("engine must be one of: auto, shapely, numpy")

    empty_xy = np.empty((0, 2), dtype=float)
    empty_i = np.empty(0, dtype=int)
    empty_f = np.empty(0, dtype=float)
    if len(i) == 0:
        return empty_xy, empty_i, empty_i.copy(), empty_i.copy(), empty_f

    proj_next = np.roll(proj, -1, axis=0)
    p = proj[i]
    r = proj_next[i] - p
    q = proj[j]
    s = proj_next[j] - q
    den = cross2(r, s)
    valid = np.abs(den) > eps
    if not np.any(valid):
        return empty_xy, empty_i, empty_i.copy(), empty_i.copy(), empty_f

    i = i[valid]
    j = j[valid]
    p = p[valid]
    r = r[valid]
    q = q[valid]
    s = s[valid]
    den = den[valid]
    qp = q - p
    t = cross2(qp, s) / den
    u = cross2(qp, r) / den
    interior = (t > eps) & (t < 1.0 - eps) & (u > eps) & (u < 1.0 - eps)
    if not np.any(interior):
        return empty_xy, empty_i, empty_i.copy(), empty_i.copy(), empty_f

    i = i[interior]
    j = j[interior]
    p = p[interior]
    r = r[interior]
    t = t[interior]
    u = u[interior]

    pts_next = np.roll(pts, -1, axis=0)
    di = pts_next[i] - pts[i]
    dj = pts_next[j] - pts[j]
    ri = pts[i] + t[:, None] * di
    rj = pts[j] + u[:, None] * dj
    # At a projected crossing, ri-rj is parallel to the viewing normal.  Thus
    # sign((di x dj).(ri-rj)) is exactly the usual oriented crossing sign.  If
    # i is the over-strand, ri-rj points toward the viewer; if j is over,
    # swapping both tangents and displacement leaves the triple product
    # unchanged.  This is also the numerator sign in the Gauss writhe
    # integral, so its spherical mean has the sign of Gauss Wr.
    triple = np.einsum("ij,ij->i", np.cross(di, dj), ri - rj)
    depth = (ri - rj) @ nvec
    nondeg = np.abs(triple) > eps

    crossing_xy = p[nondeg] + t[nondeg, None] * r[nondeg]
    signs = np.sign(triple[nondeg]).astype(int)
    return crossing_xy, signs, i[nondeg], j[nondeg], depth[nondeg]


def count_projection_crossings(
    points: np.ndarray,
    normal: Sequence[float],
    engine: str = "auto",
    eps: float = 1.0e-10,
) -> Tuple[int, int, np.ndarray]:
    """
    Count unsigned crossings and signed crossing sum for one viewing direction.

    Crossing sign is evaluated from the 3D oriented triple product
        sign[(t_i x t_j) . (r_i - r_j)]
    at the projected intersection.  This is independent of the arbitrary 2D
    screen basis, agrees with the sign of the Gauss-writhe numerator, and
    changes only at nongeneric/singular projections.

    Returns:
        unsigned_count, signed_sum, signs
    """
    _xy, signs, _i, _j, _depth = projection_crossing_details(
        points, normal, engine=engine, eps=eps
    )
    return int(len(signs)), int(signs.sum()), signs

def format_distribution(values: Iterable[int]) -> str:
    """Format integer-value frequencies in ascending order."""
    c = Counter(int(v) for v in values)
    total = sum(c.values())
    parts = []
    for key in sorted(c):
        pct = 100.0 * c[key] / total if total else 0.0
        parts.append("{}: {} ({:.2f}%)".format(key, c[key], pct))
    return "; ".join(parts)


def parse_vector(text: str) -> np.ndarray:
    """Parse 'x,y,z' into a normalized vector."""
    vals = [float(x.strip()) for x in text.split(",")]
    if len(vals) != 3:
        raise argparse.ArgumentTypeError("Direction must be x,y,z")
    return normalize_vector(vals)


def write_csv(path: str, directions: np.ndarray, unsigned: np.ndarray, signed: np.ndarray, labels: List[str]) -> None:
    """Write per-direction results."""
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["label", "nx", "ny", "nz", "unsigned_crossings", "signed_crossing_sum"])
        for label, n, nu, ws in zip(labels, directions, unsigned, signed):
            writer.writerow([label, "{:.10f}".format(n[0]), "{:.10f}".format(n[1]), "{:.10f}".format(n[2]), int(nu), int(ws)])


def run_analysis(args: argparse.Namespace) -> Tuple[str, str]:
    """Run the spherical scan and write CSV + summary files."""
    pts_full = read_xyz(args.xyz_file)
    pts = resample_closed_polyline(pts_full, args.analysis_points) if args.analysis_points else strip_closing_duplicate(pts_full)

    directions = fibonacci_sphere(args.n_directions)
    labels = ["fib_{:05d}".format(i) for i in range(len(directions))]

    if args.include_axis_views:
        axis_dirs = np.array([
            [1.0, 0.0, 0.0],  # yz projection
            [0.0, 1.0, 0.0],  # xz projection
            [0.0, 0.0, 1.0],  # xy projection
        ])
        directions = np.vstack((axis_dirs, directions))
        labels = ["yz_view(+x)", "xz_view(+y)", "xy_view(+z)"] + labels

    if args.direction:
        custom = np.vstack(args.direction)
        directions = np.vstack((custom, directions))
        labels = ["custom_{:02d}".format(i + 1) for i in range(len(custom))] + labels

    unsigned = np.zeros(len(directions), dtype=int)
    signed = np.zeros(len(directions), dtype=int)

    for k, normal in enumerate(directions):
        nu, ws, _signs = count_projection_crossings(pts, normal, engine=args.engine)
        unsigned[k] = nu
        signed[k] = ws
        if args.progress_every and ((k + 1) % args.progress_every == 0 or k + 1 == len(directions)):
            print("[INFO] analyzed {}/{} directions".format(k + 1, len(directions)))

    csv_path = args.csv_out or add_suffix(args.xyz_file, "_projection_scan", ".csv")
    summary_path = args.summary_out or add_suffix(args.xyz_file, "_projection_summary", ".txt")
    write_csv(csv_path, directions, unsigned, signed, labels)

    fib_mask = np.array([label.startswith("fib_") for label in labels], dtype=bool)
    u_scan = unsigned[fib_mask]
    s_scan = signed[fib_mask]

    min_u = int(u_scan.min())
    min_indices = np.where(fib_mask & (unsigned == min_u))[0]
    representative = min_indices[: min(args.n_representatives, len(min_indices))]

    with open(summary_path, "w", encoding="utf-8") as handle:
        handle.write("Projection-crossing analysis\n")
        handle.write("============================\n")
        handle.write("Input: {}\n".format(args.xyz_file))
        handle.write("Original points: {}\n".format(len(strip_closing_duplicate(pts_full))))
        handle.write("Analysis points: {}\n".format(len(pts)))
        handle.write("Spherical scan directions: {}\n".format(args.n_directions))
        handle.write("Crossing engine: {}\n".format(args.engine))
        handle.write("\n")
        handle.write("Mean signed crossing sum (writhe estimate): {:.8f}\n".format(float(np.mean(s_scan))))
        handle.write("Mean unsigned crossing count (ACN estimate): {:.8f}\n".format(float(np.mean(u_scan))))
        handle.write("Minimum unsigned crossing count: {}\n".format(min_u))
        handle.write("Fraction at minimum: {:.4f}%\n".format(100.0 * float(np.mean(u_scan == min_u))))
        handle.write("Signed distribution: {}\n".format(format_distribution(s_scan)))
        handle.write("Unsigned distribution: {}\n".format(format_distribution(u_scan)))
        handle.write("\nRepresentative minimum-crossing directions:\n")
        for idx in representative:
            n = directions[idx]
            handle.write(
                "  n=({:.6f}, {:.6f}, {:.6f})  unsigned={}  signed={}\n".format(
                    n[0], n[1], n[2], int(unsigned[idx]), int(signed[idx])
                )
            )
        special = [i for i, label in enumerate(labels) if not label.startswith("fib_")]
        if special:
            handle.write("\nSpecial views:\n")
            for idx in special:
                n = directions[idx]
                handle.write(
                    "  {}: n=({:.6f}, {:.6f}, {:.6f})  unsigned={}  signed={}\n".format(
                        labels[idx], n[0], n[1], n[2], int(unsigned[idx]), int(signed[idx])
                    )
                )

    print("[INFO] CSV: {}".format(csv_path))
    print("[INFO] Summary: {}".format(summary_path))
    return csv_path, summary_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Scan projection directions and count signed/unsigned crossings of a closed 3D curve.")
    parser.add_argument("xyz_file", help="Input XYZ-like coordinate file")
    parser.add_argument("--n-directions", type=int, default=10000, help="Number of Fibonacci-sphere scan directions (default: 10000)")
    parser.add_argument("--analysis-points", type=int, default=1000, help="Arc-length resample size used for the scan; 0 keeps all points (default: 1000)")
    parser.add_argument("--engine", choices=["auto", "shapely", "numpy"], default="auto", help="Crossing candidate engine (default: auto)")
    parser.add_argument("--include-axis-views", action="store_true", help="Also evaluate yz (+x), xz (+y), and xy (+z) projections")
    parser.add_argument("--direction", action="append", type=parse_vector, help="Additional viewing normal x,y,z; may be repeated")
    parser.add_argument("--n-representatives", type=int, default=8, help="Number of minimum-crossing directions to list (default: 8)")
    parser.add_argument("--progress-every", type=int, default=500, help="Print progress every N views; 0 disables (default: 500)")
    parser.add_argument("--csv-out", help="Optional output CSV path")
    parser.add_argument("--summary-out", help="Optional output summary path")
    add_gui_argument(parser)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    if gui_requested(argv):
        sys.exit(launch_gui(build_parser, __file__, argv))
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.engine == "shapely" and not HAVE_SHAPELY:
        parser.error("--engine shapely requested, but Shapely is not installed")
    run_analysis(args)


if __name__ == "__main__":
    main()
