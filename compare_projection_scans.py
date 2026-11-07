#!/usr/bin/env python3
"""
Compare direction-matched projection scans for two closed centerlines.

Example:
    python compare_projection_scans.py \
        examples/sc_W3_ALocal.xyz examples/sc_W3_ALocal_projection_scan.csv \
        examples/sc_W3_R11.xyz examples/sc_W3_R11_projection_scan.csv \
        --label-a ALocal --label-b R11 --outdir examples/conformation_comparison_figures

Run without arguments, or with --gui, to enter the parameters in a window.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np

from analyze_projection_crossings import (
    project_points,
    projection_crossing_details,
    read_xyz,
    resample_closed_polyline,
    strip_closing_duplicate,
)
from writhe_gui import add_gui_argument, gui_requested, launch_gui


CURVE_A = "#25558f"
CURVE_B = "#d16a28"
POSITIVE = "#c92b2b"
NEGATIVE = "#1768c4"
plt.rcParams["svg.fonttype"] = "none"


def safe_name(text: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in text)


def selected_formats(output_format: str) -> Tuple[str, ...]:
    return ("png", "svg") if output_format == "both" else (output_format,)


def save_figure(
    fig: plt.Figure,
    base_path: str,
    formats: Sequence[str],
    dpi: int,
) -> List[str]:
    outputs: List[str] = []
    for extension in formats:
        path = "{}.{}".format(base_path, extension)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        outputs.append(path)
    return outputs


def read_fibonacci_rows(path: str) -> Dict[str, Tuple[np.ndarray, int, int]]:
    """Read only fib_ rows as label -> (normal, unsigned, signed)."""
    rows: Dict[str, Tuple[np.ndarray, int, int]] = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            label = row["label"]
            if not label.startswith("fib_"):
                continue
            normal = np.asarray(
                [float(row["nx"]), float(row["ny"]), float(row["nz"])],
                dtype=float,
            )
            rows[label] = (
                normal,
                int(row["unsigned_crossings"]),
                int(row["signed_crossing_sum"]),
            )
    if not rows:
        raise ValueError("{} contains no fib_ rows".format(path))
    return rows


def align_rows(
    rows_a: Dict[str, Tuple[np.ndarray, int, int]],
    rows_b: Dict[str, Tuple[np.ndarray, int, int]],
) -> Dict[str, np.ndarray]:
    """Align two scans by Fibonacci label and verify their normals match."""
    labels = sorted(set(rows_a) & set(rows_b))
    if not labels:
        raise ValueError("The scans have no shared fib_ labels")
    missing_a = sorted(set(rows_b) - set(rows_a))
    missing_b = sorted(set(rows_a) - set(rows_b))
    if missing_a or missing_b:
        raise ValueError(
            "Fibonacci label sets differ: missing from A={}, missing from B={}".format(
                len(missing_a), len(missing_b)
            )
        )

    normals_a = np.vstack([rows_a[label][0] for label in labels])
    normals_b = np.vstack([rows_b[label][0] for label in labels])
    if not np.allclose(normals_a, normals_b, atol=1.0e-9, rtol=0.0):
        raise ValueError("Matched fib_ labels do not contain the same normals")

    return {
        "labels": np.asarray(labels, dtype=object),
        "normals": normals_a,
        "unsigned_a": np.asarray([rows_a[label][1] for label in labels], dtype=int),
        "signed_a": np.asarray([rows_a[label][2] for label in labels], dtype=int),
        "unsigned_b": np.asarray([rows_b[label][1] for label in labels], dtype=int),
        "signed_b": np.asarray([rows_b[label][2] for label in labels], dtype=int),
    }


def format_distribution(values: np.ndarray) -> str:
    counts = Counter(int(value) for value in values)
    total = len(values)
    return "; ".join(
        "{}: {} ({:.2f}%)".format(value, counts[value], 100.0 * counts[value] / total)
        for value in sorted(counts)
    )


def write_matched_csv(path: str, data: Dict[str, np.ndarray]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "label", "nx", "ny", "nz",
                "unsigned_a", "signed_a", "unsigned_b", "signed_b",
                "delta_unsigned_b_minus_a", "delta_signed_b_minus_a",
            ]
        )
        for index, label in enumerate(data["labels"]):
            n = data["normals"][index]
            ua = int(data["unsigned_a"][index])
            ca = int(data["signed_a"][index])
            ub = int(data["unsigned_b"][index])
            cb = int(data["signed_b"][index])
            writer.writerow(
                [label, "{:.10f}".format(n[0]), "{:.10f}".format(n[1]), "{:.10f}".format(n[2]),
                 ua, ca, ub, cb, ub - ua, cb - ca]
            )


def write_summary(
    path: str,
    data: Dict[str, np.ndarray],
    label_a: str,
    label_b: str,
) -> None:
    ua = data["unsigned_a"]
    ca = data["signed_a"]
    ub = data["unsigned_b"]
    cb = data["signed_b"]
    same_u = float(np.mean(ua == ub))
    same_c = float(np.mean(ca == cb))
    same_both = float(np.mean((ua == ub) & (ca == cb)))
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("Direction-matched conformation comparison\n")
        handle.write("========================================\n")
        handle.write("Conformation A: {}\n".format(label_a))
        handle.write("Conformation B: {}\n".format(label_b))
        handle.write("Matched Fibonacci directions: {}\n\n".format(len(ua)))
        handle.write("Metric                         {:>14} {:>14}\n".format(label_a, label_b))
        handle.write("Mean signed C(n)               {:14.6f} {:14.6f}\n".format(float(np.mean(ca)), float(np.mean(cb))))
        handle.write("Mean unsigned N(n)             {:14.6f} {:14.6f}\n".format(float(np.mean(ua)), float(np.mean(ub))))
        handle.write("Minimum unsigned N(n)          {:14d} {:14d}\n".format(int(np.min(ua)), int(np.min(ub))))
        handle.write("Fraction with N=3              {:13.2f}% {:13.2f}%\n".format(100.0 * float(np.mean(ua == 3)), 100.0 * float(np.mean(ub == 3))))
        handle.write("Fraction with C=3              {:13.2f}% {:13.2f}%\n".format(100.0 * float(np.mean(ca == 3)), 100.0 * float(np.mean(cb == 3))))
        handle.write("Fraction with C=4              {:13.2f}% {:13.2f}%\n".format(100.0 * float(np.mean(ca == 4)), 100.0 * float(np.mean(cb == 4))))
        handle.write("\nDirection-matched agreement\n")
        handle.write("Unsigned count identical: {:.2f}%\n".format(100.0 * same_u))
        handle.write("Signed count identical: {:.2f}%\n".format(100.0 * same_c))
        handle.write("Both identical: {:.2f}%\n".format(100.0 * same_both))
        handle.write("Mean |delta N|: {:.6f}\n".format(float(np.mean(np.abs(ub - ua)))))
        handle.write("Mean |delta C|: {:.6f}\n".format(float(np.mean(np.abs(cb - ca)))))
        handle.write("Maximum |delta N|: {}\n".format(int(np.max(np.abs(ub - ua)))))
        handle.write("Maximum |delta C|: {}\n".format(int(np.max(np.abs(cb - ca)))))
        handle.write("Correlation of N: {:.6f}\n".format(float(np.corrcoef(ua, ub)[0, 1])))
        handle.write("Correlation of C: {:.6f}\n".format(float(np.corrcoef(ca, cb)[0, 1])))
        handle.write("\n{} signed distribution: {}\n".format(label_a, format_distribution(ca)))
        handle.write("{} signed distribution: {}\n".format(label_b, format_distribution(cb)))
        handle.write("\n{} unsigned distribution: {}\n".format(label_a, format_distribution(ua)))
        handle.write("{} unsigned distribution: {}\n".format(label_b, format_distribution(ub)))


def distribution_plot(
    data: Dict[str, np.ndarray],
    label_a: str,
    label_b: str,
    base_path: str,
    formats: Sequence[str],
    dpi: int,
) -> List[str]:
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.4))
    for ax, key_a, key_b, title, xlabel in [
        (axes[0], "signed_a", "signed_b", "Signed crossing distribution", "C(n)"),
        (axes[1], "unsigned_a", "unsigned_b", "Unsigned crossing distribution", "N(n)"),
    ]:
        values_a = data[key_a]
        values_b = data[key_b]
        lo = int(min(values_a.min(), values_b.min()))
        hi = int(max(values_a.max(), values_b.max()))
        x = np.arange(lo, hi + 1)
        count_a = Counter(int(value) for value in values_a)
        count_b = Counter(int(value) for value in values_b)
        percent_a = np.asarray([100.0 * count_a[int(value)] / len(values_a) for value in x])
        percent_b = np.asarray([100.0 * count_b[int(value)] / len(values_b) for value in x])
        bars_a = ax.bar(x - 0.2, percent_a, width=0.4, color=CURVE_A, label=label_a)
        bars_b = ax.bar(x + 0.2, percent_b, width=0.4, color=CURVE_B, label=label_b)
        rotation = 90 if hi - lo > 8 else 0
        font_size = 5.5 if rotation else 7.0
        for bars, percentages in ((bars_a, percent_a), (bars_b, percent_b)):
            for bar, percentage in zip(bars, percentages):
                if percentage <= 0:
                    continue
                ax.annotate(
                    "{:.2f}%".format(float(percentage)),
                    xy=(bar.get_x() + bar.get_width() / 2.0, bar.get_height()),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    rotation=rotation,
                    fontsize=font_size,
                )
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("viewing directions (%)")
        tick_step = 2 if hi - lo > 12 else 1
        ax.set_xticks(np.arange(lo, hi + 1, tick_step))
        ax.set_ylim(0, max(float(percent_a.max()), float(percent_b.max())) * (1.23 if rotation else 1.16))
        ax.grid(axis="y", alpha=0.25)
        ax.legend(frameon=False)
    fig.suptitle("Same writhe, different projection-crossing distributions", fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    outputs = save_figure(fig, base_path, formats, dpi)
    plt.close(fig)
    return outputs


def difference_map_plot(
    data: Dict[str, np.ndarray],
    label_a: str,
    label_b: str,
    base_path: str,
    formats: Sequence[str],
    dpi: int,
) -> List[str]:
    normals = data["normals"]
    lon = np.arctan2(normals[:, 1], normals[:, 0])
    lat = np.arcsin(np.clip(normals[:, 2], -1.0, 1.0))
    delta_c = data["signed_b"] - data["signed_a"]
    delta_n = data["unsigned_b"] - data["unsigned_a"]
    fig = plt.figure(figsize=(12.0, 5.2))
    for index, (values, title, label) in enumerate(
        [
            (delta_c, "Signed-count difference", "delta C"),
            (delta_n, "Unsigned-count difference", "delta N"),
        ],
        start=1,
    ):
        ax = fig.add_subplot(1, 2, index, projection="mollweide")
        limit = max(1, int(np.max(np.abs(values))))
        scatter = ax.scatter(
            lon, lat, c=values, s=7, cmap="coolwarm",
            norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
        )
        ax.grid(True, alpha=0.3)
        ax.set_title("{}\n{} minus {}".format(title, label_b, label_a), fontsize=11)
        fig.colorbar(scatter, ax=ax, shrink=0.78, label=label)
    fig.suptitle("Direction-matched effect of centerline conformation", fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    outputs = save_figure(fig, base_path, formats, dpi)
    plt.close(fig)
    return outputs


def draw_projection(
    ax: plt.Axes,
    points: np.ndarray,
    normal: Sequence[float],
    title: str,
    curve_color: str,
    engine: str,
) -> None:
    projected = project_points(points, normal)
    projected = np.vstack((projected, projected[0]))
    xy, signs, _i, _j, _depth = projection_crossing_details(
        points, normal, engine=engine
    )
    ax.plot(projected[:, 0], projected[:, 1], color=curve_color, linewidth=1.7)
    if len(signs):
        colors = np.where(signs > 0, POSITIVE, NEGATIVE)
        ax.scatter(xy[:, 0], xy[:, 1], c=colors, s=32, edgecolors="white", linewidths=0.5, zorder=4)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_axis_off()
    ax.set_title("{}\nN={}, C={}".format(title, len(signs), int(signs.sum())), fontsize=11)


def coordinate_view_plot(
    xyz_a: str,
    xyz_b: str,
    label_a: str,
    label_b: str,
    analysis_points: int,
    engine: str,
    base_path: str,
    formats: Sequence[str],
    dpi: int,
) -> List[str]:
    points_a = strip_closing_duplicate(read_xyz(xyz_a))
    points_b = strip_closing_duplicate(read_xyz(xyz_b))
    if analysis_points:
        points_a = resample_closed_polyline(points_a, analysis_points)
        points_b = resample_closed_polyline(points_b, analysis_points)
    views = [("xz view", [0.0, -1.0, 0.0]), ("yz view", [1.0, 0.0, 0.0])]
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 9.0))
    for column, (view_name, normal) in enumerate(views):
        draw_projection(axes[0, column], points_a, normal, "{}: {}".format(label_a, view_name), CURVE_A, engine)
        draw_projection(axes[1, column], points_b, normal, "{}: {}".format(label_b, view_name), CURVE_B, engine)
    fig.suptitle("Coordinate-plane views of two Wr approximately +3 conformations", fontsize=15)
    fig.text(0.98, 0.015, "red = positive crossing; blue = negative crossing", ha="right", fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    outputs = save_figure(fig, base_path, formats, dpi)
    plt.close(fig)
    return outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare two direction-matched projection-crossing scans."
    )
    parser.add_argument("xyz_a", help="XYZ coordinates for conformation A")
    parser.add_argument("scan_a", help="Projection scan CSV for conformation A")
    parser.add_argument("xyz_b", help="XYZ coordinates for conformation B")
    parser.add_argument("scan_b", help="Projection scan CSV for conformation B")
    parser.add_argument("--label-a", default="conformation_A", help="Display label for A")
    parser.add_argument("--label-b", default="conformation_B", help="Display label for B")
    parser.add_argument("--analysis-points", type=int, default=1000, help="Resample size for coordinate-view drawings; 0 keeps all points")
    parser.add_argument("--engine", choices=["auto", "shapely", "numpy"], default="auto", help="Crossing engine for the coordinate-view drawings (default: auto)")
    parser.add_argument("--format", dest="output_format", choices=["png", "svg", "both"], default="png", help="Figure format: png, svg, or both (default: png)")
    parser.add_argument("--dpi", type=int, default=220, help="PNG/rasterized-element resolution (default: 220)")
    parser.add_argument("--outdir", default="conformation_comparison_figures", help="Output directory, relative to the current directory unless absolute (default: conformation_comparison_figures)")
    parser.add_argument("--stem", help="Output filename stem (default: <label-a>_vs_<label-b>)")
    add_gui_argument(parser)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    if gui_requested(argv):
        sys.exit(launch_gui(build_parser, __file__, argv))
    args = build_parser().parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    formats = selected_formats(args.output_format)
    stem = args.stem or "{}_vs_{}".format(safe_name(args.label_a), safe_name(args.label_b))
    rows_a = read_fibonacci_rows(args.scan_a)
    rows_b = read_fibonacci_rows(args.scan_b)
    data = align_rows(rows_a, rows_b)

    matched_csv = os.path.join(args.outdir, stem + "_matched_directions.csv")
    summary_path = os.path.join(args.outdir, stem + "_comparison_summary.txt")
    write_matched_csv(matched_csv, data)
    write_summary(summary_path, data, args.label_a, args.label_b)

    outputs = [matched_csv, summary_path]
    outputs.extend(
        distribution_plot(
            data, args.label_a, args.label_b,
            os.path.join(args.outdir, stem + "_distribution_comparison"),
            formats, args.dpi,
        )
    )
    outputs.extend(
        difference_map_plot(
            data, args.label_a, args.label_b,
            os.path.join(args.outdir, stem + "_direction_difference_maps"),
            formats, args.dpi,
        )
    )
    outputs.extend(
        coordinate_view_plot(
            args.xyz_a, args.xyz_b, args.label_a, args.label_b,
            args.analysis_points, args.engine,
            os.path.join(args.outdir, stem + "_coordinate_views"),
            formats, args.dpi,
        )
    )

    print("[INFO] Wrote {} output file(s) to {}".format(len(outputs), args.outdir))
    for path in outputs:
        print("  {}".format(path))


if __name__ == "__main__":
    main()
