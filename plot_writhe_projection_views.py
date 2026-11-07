#!/usr/bin/env python3
"""
Plot projection views and orientation statistics for a closed 3D curve.

Inputs:
    1. XYZ-like centerline coordinates.
    2. Optionally, the CSV produced by analyze_projection_crossings.py.
Outputs:
    PNG and/or SVG figures for the standard/representative projections and,
    when a scan CSV is supplied, histograms and spherical orientation maps.

The default representative directions reproduce the seven views discussed for
sc_W3_ALocal.xyz in the writhe report: xz, three 3-crossing oblique views, yz,
one partial-cancellation view, and one 6-crossing full-cancellation view.

Example:
    python plot_writhe_projection_views.py examples/sc_W3_ALocal.xyz \
        --scan-csv examples/sc_W3_ALocal_projection_scan.csv --annotate-crossings

Run without arguments, or with --gui, to enter the parameters in a window.

Python 3.9+; requires NumPy and Matplotlib.  Shapely 2.x is recommended because
this script imports crossing detection from analyze_projection_crossings.py.
"""

import argparse
import csv
import math
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt

from analyze_projection_crossings import (
    HAVE_SHAPELY,
    count_projection_crossings,
    project_points,
    projection_crossing_details,
    read_xyz,
    resample_closed_polyline,
    strip_closing_duplicate,
)

try:
    from writhe_gui import add_gui_argument, gui_requested, launch_gui
except ModuleNotFoundError as exc:  # the command line works without the GUI module
    if exc.name != "writhe_gui":
        raise

    def add_gui_argument(parser: argparse.ArgumentParser) -> None:
        """writhe_gui.py is missing, so there is no --gui option."""

    def gui_requested(argv: Sequence[str]) -> bool:
        return False

    launch_gui = None


DEFAULT_VIEWS = [
    ("xz_reference", np.array([0.0, -1.0, 0.0])),
    ("representative_A", np.array([0.2494, 0.8564, 0.4521])),
    ("representative_B", np.array([-0.2700, 0.8643, 0.4243])),
    ("representative_C", np.array([-0.0244, 0.6586, 0.7521])),
    ("yz_reference", np.array([1.0, 0.0, 0.0])),
    ("partial_cancellation", np.array([-0.1754, -0.6293, 0.7571])),
    ("full_cancellation", np.array([0.1941, 0.2748, 0.9417])),
]

CURVE_COLOR = "#25558f"
POSITIVE_COLOR = "#c92b2b"
NEGATIVE_COLOR = "#1768c4"

# Preserve labels as editable SVG text instead of converting every glyph to a
# vector path. Viewers will use the declared sans-serif font or a substitute.
plt.rcParams["svg.fonttype"] = "none"


def normalize(v: Sequence[float]) -> np.ndarray:
    arr = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(arr))
    if n <= 1.0e-15:
        raise ValueError("View normal cannot be zero.")
    return arr / n


def parse_named_view(text: str) -> Tuple[str, np.ndarray]:
    """Parse NAME:x,y,z from the command line."""
    if ":" not in text:
        raise argparse.ArgumentTypeError("Use NAME:x,y,z")
    name, vec = text.split(":", 1)
    vals = [float(x.strip()) for x in vec.split(",")]
    if len(vals) != 3:
        raise argparse.ArgumentTypeError("Use NAME:x,y,z")
    return name.strip(), normalize(vals)


def safe_name(text: str) -> str:
    """Convert an arbitrary label into a filename-safe stem."""
    return "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in text)


def output_prefix(xyz_file: str, outdir: Optional[str]) -> Tuple[str, str]:
    """Return output directory and input basename without extension."""
    root = os.path.splitext(os.path.basename(xyz_file))[0]
    directory = outdir or os.path.join(os.path.dirname(os.path.abspath(xyz_file)), root + "_writhe_figures")
    os.makedirs(directory, exist_ok=True)
    return directory, root


def selected_formats(output_format: str) -> Tuple[str, ...]:
    """Expand the CLI output-format selection into filename extensions."""
    if output_format == "both":
        return ("png", "svg")
    return (output_format,)


def save_figure(
    fig: plt.Figure,
    base_path: str,
    formats: Sequence[str],
    dpi: int,
) -> List[str]:
    """Save one Matplotlib figure in every requested format."""
    outputs: List[str] = []
    for extension in formats:
        path = "{}.{}".format(base_path, extension)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        outputs.append(path)
    return outputs


def plot_one_projection(
    points: np.ndarray,
    normal: Sequence[float],
    title: str,
    base_path: str,
    engine: str,
    annotate_crossings: bool,
    dpi: int,
    formats: Sequence[str],
) -> List[str]:
    """Plot one orthographic view and optionally mark crossing signs."""
    n = normalize(normal)
    proj = project_points(points, n)
    closed_proj = np.vstack((proj, proj[0]))
    crossing_xy, signs, _i, _j, _depth = projection_crossing_details(points, n, engine=engine)
    unsigned = int(len(signs))
    signed = int(signs.sum())

    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    ax.plot(closed_proj[:, 0], closed_proj[:, 1], linewidth=2.0, color=CURVE_COLOR)
    if annotate_crossings and unsigned:
        colors = np.where(signs > 0, POSITIVE_COLOR, NEGATIVE_COLOR)
        ax.scatter(crossing_xy[:, 0], crossing_xy[:, 1], s=52, c=colors, edgecolors="white", linewidths=0.6, zorder=5)
        for xy, sign in zip(crossing_xy, signs):
            ax.text(xy[0], xy[1], "+" if sign > 0 else "-", color="white", ha="center", va="center", fontsize=7, zorder=6)

    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("projection coordinate 1")
    ax.set_ylabel("projection coordinate 2")
    ax.set_title(
        "{}\nnormal = ({:.4f}, {:.4f}, {:.4f}); crossings = {}; signed sum = {}".format(
            title, n[0], n[1], n[2], unsigned, signed
        )
    )
    fig.tight_layout()
    outputs = save_figure(fig, base_path, formats, dpi)
    plt.close(fig)
    return outputs


def draw_projection_panel(
    ax: plt.Axes,
    points: np.ndarray,
    name: str,
    normal: Sequence[float],
    engine: str,
) -> None:
    """Draw one compact, sign-annotated projection in a multipanel figure."""
    n = normalize(normal)
    proj = project_points(points, n)
    closed_proj = np.vstack((proj, proj[0]))
    crossing_xy, signs, _i, _j, _depth = projection_crossing_details(points, n, engine=engine)
    ax.plot(closed_proj[:, 0], closed_proj[:, 1], linewidth=1.5, color=CURVE_COLOR)
    if len(signs):
        colors = np.where(signs > 0, POSITIVE_COLOR, NEGATIVE_COLOR)
        ax.scatter(crossing_xy[:, 0], crossing_xy[:, 1], s=31, c=colors, edgecolors="white", linewidths=0.45, zorder=5)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_axis_off()
    label = name.replace("_", " ")
    ax.set_title(
        "{}\n$n$=({:.4f}, {:.4f}, {:.4f})\n$N$={}, $C$={}".format(
            label, n[0], n[1], n[2], len(signs), int(signs.sum())
        ),
        fontsize=9,
    )


def plot_report_panels(
    points: np.ndarray,
    views: Sequence[Tuple[str, np.ndarray]],
    outdir: str,
    stem: str,
    engine: str,
    dpi: int,
    formats: Sequence[str],
) -> List[str]:
    """Generate the two coordinate-view and seven-view report panels."""
    by_name = {name: normal for name, normal in views}
    required = [name for name, _normal in DEFAULT_VIEWS]
    missing = [name for name in required if name not in by_name]
    if missing:
        raise ValueError("--report-panels requires the built-in views: {}".format(", ".join(missing)))

    outputs: List[str] = []
    coordinate_names = ["xz_reference", "yz_reference"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 5.6))
    for ax, name in zip(axes, coordinate_names):
        draw_projection_panel(ax, points, name, by_name[name], engine)
    fig.suptitle("The same centerline has different projected crossing counts", fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    base_path = os.path.join(outdir, "{}_coordinate_plane_comparison".format(stem))
    outputs.extend(save_figure(fig, base_path, formats, dpi))
    plt.close(fig)

    fig, axes = plt.subplots(2, 4, figsize=(14.5, 8.2))
    for ax, name in zip(axes.flat, required):
        draw_projection_panel(ax, points, name, by_name[name], engine)
    axes.flat[-1].set_axis_off()
    fig.suptitle("Seven representative projections of the same positive-writhe centerline", fontsize=16)
    fig.text(
        0.99,
        0.015,
        "red = positive crossing    blue = negative crossing",
        ha="right",
        fontsize=10,
        color="#333333",
    )
    fig.tight_layout(rect=(0, 0.035, 1, 0.94))
    base_path = os.path.join(outdir, "{}_representative_projections".format(stem))
    outputs.extend(save_figure(fig, base_path, formats, dpi))
    plt.close(fig)
    return outputs


def read_scan_csv(path: str) -> Dict[str, np.ndarray]:
    """Read the CSV produced by analyze_projection_crossings.py."""
    labels: List[str] = []
    vecs: List[List[float]] = []
    unsigned: List[int] = []
    signed: List[int] = []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            labels.append(row["label"])
            vecs.append([float(row["nx"]), float(row["ny"]), float(row["nz"])])
            unsigned.append(int(row["unsigned_crossings"]))
            signed.append(int(row["signed_crossing_sum"]))
    return {
        "labels": np.asarray(labels, dtype=object),
        "directions": np.asarray(vecs, dtype=float),
        "unsigned": np.asarray(unsigned, dtype=int),
        "signed": np.asarray(signed, dtype=int),
    }


def scan_only(data: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Remove axis/custom rows and keep Fibonacci-scan rows when available."""
    labels = data["labels"]
    mask = np.array([str(x).startswith("fib_") for x in labels], dtype=bool)
    if not np.any(mask):
        mask = np.ones(len(labels), dtype=bool)
    return {key: value[mask] for key, value in data.items()}


def plot_histograms(
    data: Dict[str, np.ndarray],
    outdir: str,
    stem: str,
    dpi: int,
    formats: Sequence[str],
) -> List[str]:
    """Plot signed and unsigned crossing distributions."""
    data = scan_only(data)
    outputs: List[str] = []
    for field, xlabel, suffix in [
        ("signed", "signed crossing sum C(n)", "signed_histogram"),
        ("unsigned", "unsigned crossing count N(n)", "unsigned_histogram"),
    ]:
        values = data[field]
        categories, counts = np.unique(values, return_counts=True)
        percentages = 100.0 * counts / len(values)
        figure_width = min(11.0, max(7.0, 3.8 + 0.38 * len(categories)))
        fig, ax = plt.subplots(figsize=(figure_width, 4.8))
        bars = ax.bar(categories, counts, width=0.82)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("number of viewing directions")
        ax.set_title("{}: mean = {:.4f}".format(xlabel, float(np.mean(values))))
        rotation = 90 if len(categories) > 8 else 0
        font_size = 6.5 if len(categories) > 12 else 8.0
        for bar, percentage in zip(bars, percentages):
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
        tick_values = categories
        if len(categories) > 14:
            tick_values = categories[::2]
            if tick_values[-1] != categories[-1]:
                tick_values = np.append(tick_values, categories[-1])
        ax.set_xticks(tick_values)
        ax.set_ylim(0, max(counts) * (1.23 if rotation else 1.16))
        ax.grid(axis="y", alpha=0.2)
        fig.tight_layout()
        base_path = os.path.join(outdir, "{}_{}".format(stem, suffix))
        outputs.extend(save_figure(fig, base_path, formats, dpi))
        plt.close(fig)
    return outputs


def plot_orientation_maps(
    data: Dict[str, np.ndarray],
    outdir: str,
    stem: str,
    dpi: int,
    formats: Sequence[str],
) -> List[str]:
    """Plot Mollweide maps of directional signed and unsigned crossings."""
    data = scan_only(data)
    n = data["directions"]
    lon = np.arctan2(n[:, 1], n[:, 0])
    lat = np.arcsin(np.clip(n[:, 2], -1.0, 1.0))
    outputs: List[str] = []
    for field, title, suffix in [
        ("signed", "Directional algebraic crossing sum C(n)", "signed_orientation_map"),
        ("unsigned", "Unsigned crossing count N(n)", "unsigned_orientation_map"),
    ]:
        fig = plt.figure(figsize=(9.0, 5.2))
        ax = fig.add_subplot(111, projection="mollweide")
        color_args = {}
        if field == "unsigned":
            # Rare near-polar projections have much larger counts.  Clipping
            # the display at 10 preserves the report's readable color scale;
            # the CSV and histogram retain the full values.
            color_args = {"vmin": max(0, int(np.min(data[field]))), "vmax": 10}
        sc = ax.scatter(lon, lat, c=data[field], s=8, **color_args)
        ax.grid(True, alpha=0.35)
        ax.set_title(title)
        fig.colorbar(sc, ax=ax, label=field)
        fig.tight_layout()
        base_path = os.path.join(outdir, "{}_{}".format(stem, suffix))
        outputs.extend(save_figure(fig, base_path, formats, dpi))
        plt.close(fig)
    return outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot representative projection views and scan statistics for a closed curve.")
    parser.add_argument("xyz_file", help="Input XYZ-like coordinate file")
    parser.add_argument("--scan-csv", help="CSV from analyze_projection_crossings.py; enables histograms and orientation maps")
    parser.add_argument("--analysis-points", type=int, default=0, help="Resample to this many points for plotting/crossing annotation; 0 uses all points")
    parser.add_argument("--engine", choices=["auto", "shapely", "numpy"], default="auto", help="Crossing engine (default: auto)")
    parser.add_argument("--annotate-crossings", action="store_true", help="Mark each crossing with + or -")
    parser.add_argument("--view", action="append", type=parse_named_view, help="Add a custom NAME:x,y,z view; may be repeated")
    parser.add_argument("--only-custom-views", action="store_true", help="Do not include the seven built-in report views")
    parser.add_argument("--report-panels", action="store_true", help="Also generate combined coordinate-plane and seven-view report panels")
    parser.add_argument("--outdir", help="Output directory; default is <input>_writhe_figures")
    parser.add_argument("--format", dest="output_format", choices=["png", "svg", "both"], default="png", help="Figure format: png, svg, or both (default: png)")
    parser.add_argument("--dpi", type=int, default=220, help="PNG/rasterized-element resolution (default: 220); SVG geometry remains vector")
    add_gui_argument(parser)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    if gui_requested(argv):
        sys.exit(launch_gui(build_parser, __file__, argv))
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "gui", False):  # an abbreviation of --gui, such as --gu
        sys.exit(launch_gui(build_parser, __file__, argv))
    if args.engine == "shapely" and not HAVE_SHAPELY:
        parser.error("--engine shapely requested, but Shapely is not installed")
    pts = strip_closing_duplicate(read_xyz(args.xyz_file))
    if args.analysis_points:
        pts = resample_closed_polyline(pts, args.analysis_points)

    outdir, stem = output_prefix(args.xyz_file, args.outdir)
    formats = selected_formats(args.output_format)
    views: List[Tuple[str, np.ndarray]] = [] if args.only_custom_views else list(DEFAULT_VIEWS)
    if args.view:
        views.extend(args.view)

    outputs: List[str] = []
    for name, normal in views:
        base_path = os.path.join(outdir, "{}_{}".format(stem, safe_name(name)))
        outputs.extend(
            plot_one_projection(
                pts,
                normal,
                name,
                base_path,
                engine=args.engine,
                annotate_crossings=args.annotate_crossings,
                dpi=args.dpi,
                formats=formats,
            )
        )

    if args.report_panels:
        outputs.extend(
            plot_report_panels(
                pts, views, outdir, stem, args.engine, args.dpi, formats
            )
        )

    if args.scan_csv:
        data = read_scan_csv(args.scan_csv)
        outputs.extend(plot_histograms(data, outdir, stem, args.dpi, formats))
        outputs.extend(plot_orientation_maps(data, outdir, stem, args.dpi, formats))

    print("[INFO] Wrote {} figure(s) to {}".format(len(outputs), outdir))
    for path in outputs:
        print("  {}".format(path))


if __name__ == "__main__":
    main()
