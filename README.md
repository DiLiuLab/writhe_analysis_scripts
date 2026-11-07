# Writhe Projection Analysis Scripts

These scripts reproduce the numerical analysis used for the
`sc_W3_ALocal.xyz` and `sc_W3_R11.xyz` writhe/projection examples. They
distinguish the **directional signed crossing sum** `w(n)` from the **unsigned
crossing count** `N(n)`, estimate 3D writhe by averaging `w(n)` over many
viewing directions, and compare two conformations on exactly the same sampled
viewing normals.

Throughout this README, `C(n)` and `w(n)` mean the same directional signed
crossing sum. The reports use `C(n)` to keep it visually distinct from the
three-dimensional writhe `Wr`.

Every script can be used from the command line or through a small graphical
interface (see [Graphical interface](#graphical-interface)).

## Repository layout

| Path | Contents |
|---|---|
| `analyze_projection_crossings.py` | Stage 1: scan viewing directions and count signed/unsigned crossings |
| `plot_writhe_projection_views.py` | Stage 2: projection figures, histograms, and orientation maps |
| `validate_projection_scan.py` | Stage 3: recheck a reduced-resolution scan at full resolution |
| `compare_projection_scans.py` | Stage 4: direction-matched comparison of two conformations |
| `writhe_gui.py` | Graphical interface shared by the four scripts; run it directly for one window with all four |
| `test_crossing_signs.py`, `test_gui.py` | Unit tests |
| `requirements.txt` | Required Python packages |
| `examples/` | Example input centerlines and every output used in the report |

The `examples/` directory contains:

| File or folder | Contents |
|---|---|
| `sc_W3_ALocal.xyz`, `sc_W3_R11.xyz` | The two report centerlines: 2,000 points, closed length 1,071 Å, `Wr` approximately `+3` |
| `sc_L714_Wm2_Logo.xyz` | A negative-writhe example: 2,000 points, closed length 714 Å, `Wr` approximately `-2` (see [Negative-writhe example](#negative-writhe-example-sc_l714_wm2_logoxyz)) |
| `sc_W3_ALocal_fake.pdb` | `sc_W3_ALocal.xyz` written as a pseudo-PDB file (one atom per point) for molecular viewers |
| `sc_xyz_generation_parameters.md` | How the two `sc_W3_*` centerlines were generated |
| `*_projection_scan.csv`, `*_projection_summary.txt` | 10,000-direction scans at 1,000 analysis points, one per centerline |
| `*_projection_scan_validation.txt` | Full-resolution validation reports for those scans |
| `sc_W3_ALocal_writhe_figures/`, `sc_W3_R11_writhe_figures/`, `sc_L714_Wm2_Logo_writhe_figures/` | Figures from the plotting script |
| `conformation_comparison_figures/` | Direction-matched comparison of the two `sc_W3_*` conformations |

## Installation

Python 3.9 or later is sufficient. Clone the repository and install the
required packages:

```bash
git clone https://github.com/DiLiuLab/writhe_analysis_scripts.git
cd writhe_analysis_scripts
pip install -r requirements.txt
```

`requirements.txt` lists NumPy, Matplotlib, and Shapely. `Shapely >= 2` is
strongly recommended because its spatial index makes the crossing search much
faster. `analyze_projection_crossings.py` also has a NumPy-only fallback
(`--engine numpy`), but it is much slower for large curves or many
orientations.

Confirm the installed packages with:

```bash
python -c "import numpy, matplotlib; print(numpy.__version__, matplotlib.__version__)"
python -c "import shapely; print(shapely.__version__)"
```

If the Shapely import fails, either install it or explicitly use
`--engine numpy`. On systems where Matplotlib cannot write its normal cache
directory, a temporary writable cache may be selected for the command:

```bash
MPLCONFIGDIR=/tmp/writhe_matplotlib \
python plot_writhe_projection_views.py examples/sc_W3_ALocal.xyz \
    --scan-csv examples/sc_W3_ALocal_projection_scan.csv
```

The graphical interface needs only `tkinter`, which is part of the Python
standard library. The python.org and Anaconda installers include it; with
Homebrew Python install `python-tk`, and on Debian/Ubuntu install
`python3-tk`.

## Quick start

Run the commands in this README from the repository root (the directory that
contains the scripts). The example inputs are in `examples/`.

The normal workflow has four stages:

| Stage | Script | Main input | Main output |
|---|---|---|---|
| 1. Analyze | `analyze_projection_crossings.py` | XYZ centerline | Scan CSV and text summary |
| 2. Draw | `plot_writhe_projection_views.py` | XYZ and, optionally, scan CSV | PNG and/or SVG projections, histograms, and orientation maps |
| 3. Validate | `validate_projection_scan.py` | XYZ and scan CSV | Full-resolution comparison report |
| 4. Compare | `compare_projection_scans.py` | Two XYZ files and their scan CSVs | Matched-direction tables and PNG/SVG comparison figures |

Minimal report workflow:

```bash
python analyze_projection_crossings.py examples/sc_W3_ALocal.xyz \
    --n-directions 10000 \
    --analysis-points 1000

python plot_writhe_projection_views.py examples/sc_W3_ALocal.xyz \
    --scan-csv examples/sc_W3_ALocal_projection_scan.csv \
    --annotate-crossings \
    --report-panels

python validate_projection_scan.py examples/sc_W3_ALocal.xyz \
    examples/sc_W3_ALocal_projection_scan.csv \
    --n-check 500 \
    --seed 2026

python analyze_projection_crossings.py examples/sc_W3_R11.xyz \
    --n-directions 10000 \
    --analysis-points 1000

python compare_projection_scans.py \
    examples/sc_W3_ALocal.xyz examples/sc_W3_ALocal_projection_scan.csv \
    examples/sc_W3_R11.xyz examples/sc_W3_R11_projection_scan.csv \
    --label-a "ALocal" \
    --label-b "R11" \
    --format both \
    --outdir examples/conformation_comparison_figures
```

The backslash (`\`) at the end of a shell-command line means “continue this
command on the next line.” It is not part of an option value. The same command
can be written on one long line without backslashes.

These commands write their results into `examples/`, replacing the copies
that come with the repository. To get the original files back, run
`git checkout -- examples`.

## Graphical interface

Each script opens a window for its parameters when it is started with no
arguments or with `--gui`:

```bash
python analyze_projection_crossings.py
python plot_writhe_projection_views.py --gui
```

Arguments given together with `--gui` pre-fill the form, which is a convenient
way to start from one of the commands in this README and adjust it:

```bash
python analyze_projection_crossings.py examples/sc_W3_ALocal.xyz --n-directions 1000 --gui
```

To get all four scripts as tabs of one window, in workflow order (Analyze,
Plot, Validate, Compare), run:

```bash
python writhe_gui.py
```

The window lists every command-line option of the script with its `--help`
text and default value. Required fields are marked `*`, and file and folder
fields have a **Browse…** button. As you edit, the **Equivalent command** box
shows the command line that does the same thing; **Copy command** puts it on
the clipboard so a GUI run can be repeated in a terminal or a batch script.
**Run** checks the values with the script's own command-line parser, then
runs the script as a separate Python process and shows its messages in the
log. **Stop** ends a long scan.

- Relative paths are resolved from the working folder shown at the top of the
  window, which is the directory the GUI was started from. A file chosen with
  **Browse…** inside that folder is entered as a relative path.
- Values left at their defaults are omitted from the command, as you would
  omit them when typing it. Clearing a field restores its default.
- Repeatable options such as `--direction` and `--view` take one value per
  line. Vectors beginning with a minus sign need no special care; the GUI
  writes them as `--direction=-0.27,0.86,0.42`.
- Command-line use is unchanged: whenever arguments are given without
  `--gui`, the script runs directly. If no window can be opened (for example,
  over SSH or when `tkinter` is missing), a script started without arguments
  prints the reason and its `--help` text instead.

## Input coordinate format

The XYZ reader uses the first three numeric fields from every nonblank,
noncomment line. Lines beginning with `#` or `!` are ignored. Extra columns are
allowed. If the last coordinate duplicates the first coordinate, the duplicate
is removed; the scripts then treat the polyline as closed by connecting the
last retained point back to the first.

Example:

```text
# x y z
12.0  0.0  1.5
11.8  0.7  1.9
...
```

At least four valid three-dimensional points are required. Coordinates must be
finite numbers.

## Files

### `analyze_projection_crossings.py`

**Function.** Reads a closed 3D centerline, scans viewing directions on a Fibonacci sphere, detects regular self-crossings in each orthographic projection, assigns a crossing sign, and writes per-direction signed and unsigned crossing counts.

Important functions:

- `read_xyz()` — reads an XYZ-like coordinate file.
- `resample_closed_polyline()` — resamples the curve uniformly in arc length for faster scans.
- `fibonacci_sphere()` — generates nearly uniform viewing directions.
- `project_points()` — orthographically projects the 3D curve onto a plane normal to a chosen viewing direction.
- `projection_crossing_details()` — finds all regular crossings and their signs for one projection.
- `count_projection_crossings()` — returns the unsigned crossing count and signed crossing sum for one projection.
- `run_analysis()` — performs the full orientation scan and writes the CSV and text summary.

**Typical use for the report example:**

```bash
python analyze_projection_crossings.py examples/sc_W3_ALocal.xyz \
    --n-directions 10000 \
    --analysis-points 1000 \
    --include-axis-views \
    --direction 0.2494,0.8564,0.4521 \
    --direction=-0.2700,0.8643,0.4243 \
    --direction=-0.0244,0.6586,0.7521 \
    --direction=-0.1754,-0.6293,0.7571 \
    --direction 0.1941,0.2748,0.9417
```

The negative-vector options use `--direction=...` because a value beginning with `-` can otherwise be interpreted as a command-line option.

#### Analysis options

| Option | Meaning |
|---|---|
| `xyz_file` | Required XYZ-like coordinate file. |
| `--n-directions N` | Number of nearly uniform Fibonacci-sphere normals. Default: `10000`. Only these directions receive `fib_` labels and contribute to the reported spherical statistics. |
| `--analysis-points N` | Resample the closed centerline to `N` equally spaced arc-length points before scanning. Default: `1000`. Use `0` to keep all original points. |
| `--engine auto` | Use Shapely when available and otherwise use the NumPy fallback. |
| `--engine shapely` | Require the faster Shapely spatial-index implementation. |
| `--engine numpy` | Force the portable NumPy implementation. It can be much slower. |
| `--include-axis-views` | Add the yz (`+x` normal), xz (`+y` normal), and xy (`+z` normal) projections to the CSV as diagnostic rows. These are not included in the spherical average. |
| `--direction x,y,z` | Add one explicitly requested normal. Repeat the option for multiple normals. The vector is normalized internally. |
| `--n-representatives N` | List up to `N` minimum-unsigned-crossing Fibonacci directions in the summary. Default: `8`. |
| `--progress-every N` | Print progress after every `N` directions. Use `0` to disable. |
| `--csv-out PATH` | Override the default scan-CSV filename. |
| `--summary-out PATH` | Override the default text-summary filename. |
| `--gui` | Open the graphical interface, pre-filled with any other arguments given. |

#### What `--analysis-points` changes

For `sc_W3_ALocal.xyz`, the original centerline contains 2,000 points. The
option

```bash
--analysis-points 1000
```

constructs a new 1,000-point closed polyline with approximately equal
arc-length spacing and performs all scan crossing calculations on that
polyline. This substantially reduces the cost of a 10,000-direction scan.
Using `--analysis-points 0` retains all 2,000 original points and is slower but
can resolve a direction lying very close to a crossing-change boundary more
faithfully. Resampling does not edit the XYZ file.

#### CSV labels and columns

The output CSV contains:

| Column | Meaning |
|---|---|
| `label` | Origin of the direction: `fib_#####`, `custom_##`, or an axis-view label. |
| `nx`, `ny`, `nz` | Components of the normalized viewing direction. |
| `unsigned_crossings` | Total visible proper crossings `N(n)`. |
| `signed_crossing_sum` | Algebraic crossing sum `C(n)`. |

Rows labeled `fib_00000`, `fib_00001`, and so on are the uniform spherical
sample. `custom_` and axis-view rows are useful diagnostics but do not
contribute to the mean, distribution, minimum, or fraction reported in the
summary file.

Outputs are written next to the input file, with names derived from the input
name:

- `examples/sc_W3_ALocal_projection_scan.csv`
- `examples/sc_W3_ALocal_projection_summary.txt`

For this curve, the report-level result is approximately:

- spherical mean signed crossing sum: `+3.0` (the writhe estimate)
- spherical mean unsigned crossing count: `4.88` (an ACN estimate, not writhe)
- smallest unsigned count found: `3`
- xz projection: 3 crossings, signed sum +3
- yz projection: 4 crossings, signed sum +4
- `n=(0.1941,0.2748,0.9417)`: 6 crossings, signed sum 0

A quick 1000-direction run with the packaged script gives mean signed `+3.002`, mean unsigned `4.882`, minimum `3`, and 14.7% of directions at the minimum.

For a faster exploratory scan, use fewer directions and fewer curve points:

```bash
python analyze_projection_crossings.py examples/sc_W3_ALocal.xyz \
    --n-directions 1000 \
    --analysis-points 500
```

Use the larger settings for the report-quality calculation.

### `plot_writhe_projection_views.py`

**Function.** Generates the projection figures used to explain the problem visually. By default it plots the xz and yz views, the three representative 3-crossing oblique orientations, and the partial- and full-cancellation examples. If supplied with the scan CSV, it also generates signed/unsigned histograms and spherical orientation maps. `--report-panels` additionally writes the two-panel coordinate comparison and the seven-panel representative-projection figure used in the reports. Output can be PNG, SVG, or both.

Important functions:

- `plot_one_projection()` — draws one orthographic projection and optionally labels each crossing `+` or `-`.
- `plot_histograms()` — plots the directional distributions of signed and unsigned crossings.
- `plot_orientation_maps()` — maps crossing counts versus viewing direction on a Mollweide projection of the viewing sphere.

**Typical use:**

```bash
python plot_writhe_projection_views.py examples/sc_W3_ALocal.xyz \
    --scan-csv examples/sc_W3_ALocal_projection_scan.csv \
    --analysis-points 1000 \
    --annotate-crossings \
    --report-panels
```

Unless `--outdir` is supplied, the figures are written to a directory named
after the input, next to the input file: here,
`examples/sc_W3_ALocal_writhe_figures/`.

#### Plotting options and what they affect

| Option | Meaning |
|---|---|
| `xyz_file` | Required even when plotting CSV statistics, because the script can also draw projection views. |
| `--scan-csv PATH` | Read the scan data and generate the signed/unsigned histograms and orientation maps. Without this option, no histogram or sphere map is made. |
| `--analysis-points N` | Resolution used only for newly calculated projection drawings and crossing annotations. It does not recompute or change values already stored in the CSV. Default: `0`, meaning all XYZ points. |
| `--annotate-crossings` | Mark crossings in individual projection images with red positive or blue negative symbols. It does not alter histograms or maps. |
| `--view NAME:x,y,z` | Add a named custom projection. Repeat for more views. For a vector beginning with a minus sign, use `--view=NAME:-x,y,z`. |
| `--only-custom-views` | Suppress the seven built-in individual projection images. This is useful when only CSV-derived histograms/maps or explicitly supplied custom views are wanted. |
| `--report-panels` | Additionally create the combined xz/yz panel and combined seven-view panel used in the reports. Do not combine this with `--only-custom-views`, because the report panels require the built-in views. |
| `--outdir PATH` | Choose the output directory. |
| `--format png` | Write PNG files. This is the default and preserves the earlier script behavior. |
| `--format svg` | Write scalable vector graphics instead of PNG files. |
| `--format both` | Write both `.png` and `.svg` versions of every requested figure. |
| `--dpi N` | Resolution for PNG output and any rasterized SVG elements. Default: `220`. Lines, markers, and text geometry in SVG remain vector-based. |
| `--gui` | Open the graphical interface, pre-filled with any other arguments given. |

#### SVG output

Generate the standard views, histograms, and maps as SVG:

```bash
python plot_writhe_projection_views.py examples/sc_W3_ALocal.xyz \
    --scan-csv examples/sc_W3_ALocal_projection_scan.csv \
    --format svg \
    --outdir examples/sc_W3_ALocal_writhe_figures
```

Generate both PNG and SVG versions:

```bash
python plot_writhe_projection_views.py examples/sc_W3_ALocal.xyz \
    --scan-csv examples/sc_W3_ALocal_projection_scan.csv \
    --annotate-crossings \
    --report-panels \
    --format both \
    --outdir examples/sc_W3_ALocal_writhe_figures
```

SVG is useful for manuscripts, presentations, and later editing in programs
such as Inkscape or Adobe Illustrator. Projection curves, text, axes, and
markers remain scalable, and labels are retained as editable SVG text rather
than being converted into letter-shaped paths. The 10,000-point orientation
maps are also saved as vector scatter plots, so their SVG files are
substantially larger than the histogram SVG files. For this supplied scan,
each orientation-map SVG is approximately a few megabytes. A viewer that does
not have the declared sans-serif font may substitute a similar font.

#### Which CSV rows are plotted?

The histogram and orientation-map functions first look for labels beginning
with `fib_`. If at least one exists, they use **only** the `fib_` rows. Axis and
custom rows are excluded so they do not bias the viewing-sphere distribution.
If a user-supplied CSV contains no `fib_` labels, the plotting script falls
back to using all rows.

The projection-view drawings do not use CSV labels. They are calculated
directly from the XYZ coordinates and the specified built-in or custom viewing
normals.

#### Generate histograms and orientation maps without the standard views

```bash
python plot_writhe_projection_views.py examples/sc_W3_ALocal.xyz \
    --scan-csv examples/sc_W3_ALocal_projection_scan.csv \
    --only-custom-views \
    --outdir examples/sc_W3_ALocal_writhe_figures
```

`--only-custom-views` is optional. It merely avoids regenerating the seven
standard individual projections. The current program creates the two
histograms and the two orientation maps together whenever `--scan-csv` is
provided; it does not currently have a histogram-only flag.

Every nonzero histogram bar is labeled with its percentage of the Fibonacci
viewing-direction sample. The vertical axis remains the absolute number of
directions, so both counts and percentages can be read from the same figure.

The CSV-derived files use the selected extension (`.png`, `.svg`, or both):

- `sc_W3_ALocal_signed_histogram.<format>`
- `sc_W3_ALocal_unsigned_histogram.<format>`
- `sc_W3_ALocal_signed_orientation_map.<format>`
- `sc_W3_ALocal_unsigned_orientation_map.<format>`

#### Generate the report projection panels as well

```bash
python plot_writhe_projection_views.py examples/sc_W3_ALocal.xyz \
    --scan-csv examples/sc_W3_ALocal_projection_scan.csv \
    --analysis-points 1000 \
    --annotate-crossings \
    --report-panels \
    --outdir examples/sc_W3_ALocal_writhe_figures
```

This generates the seven individual projection PNGs, the two combined report
panels, both histograms, and both orientation maps.

#### Report figure-to-file mapping

The ALocal figures are in `examples/sc_W3_ALocal_writhe_figures/` and the R11
figures in `examples/sc_W3_R11_writhe_figures/`.

| Report figure | Generated file |
|---|---|
| Coordinate-plane comparison (Figure 4) | `sc_W3_ALocal_coordinate_plane_comparison.<format>` |
| Seven representative views (Figure 5) | `sc_W3_ALocal_representative_projections.<format>` |
| Unsigned orientation map (Figure 6) | `sc_W3_ALocal_unsigned_orientation_map.<format>` |
| Unsigned histogram (Figure 7) | `sc_W3_ALocal_unsigned_histogram.<format>` |
| Signed orientation map (Figure 8A) | `sc_W3_ALocal_signed_orientation_map.<format>` |
| Signed C(n) histogram (Figure 8B) | `sc_W3_ALocal_signed_histogram.<format>` |
| R11 unsigned orientation map (Figure 12) | `sc_W3_R11_unsigned_orientation_map.<format>` |
| R11 signed orientation map (Figure 13) | `sc_W3_R11_signed_orientation_map.<format>` |
| R11 signed C(n) histogram (Figure 14) | `sc_W3_R11_signed_histogram.<format>` |

Each point in an orientation map represents one `fib_` viewing normal. The
longitude and latitude are computed as `atan2(ny,nx)` and `asin(nz)`, and the
point color is the corresponding `N(n)` or `C(n)` value from the CSV. The map
is a scatter plot of the sampled directions; it does not interpolate a
continuous surface between samples. The unsigned map clips its displayed
color scale at `10` for readability, while the CSV and histogram retain the
full unsigned counts.

You can add another projection with:

```bash
python plot_writhe_projection_views.py examples/sc_W3_ALocal.xyz \
    --view my_view:0.2,0.3,0.9
```

To draw only that custom view, add `--only-custom-views`:

```bash
python plot_writhe_projection_views.py examples/sc_W3_ALocal.xyz \
    --only-custom-views \
    --view my_view:0.2,0.3,0.9 \
    --annotate-crossings
```

### `validate_projection_scan.py`

**Function.** Checks that a reduced-resolution spherical scan is not creating or losing crossings. It randomly selects directions from the scan CSV and recomputes them using the original full-resolution centerline. It can also run a separate, larger coarse orientation search to look for a missed low-crossing region.

Important functions:

- `validate_rows()` — recomputes selected viewing directions at full coordinate resolution and reports mismatches.
- `coarse_minimum_search()` — performs an independent spherical search for the minimum unsigned crossing count.

**Full-resolution validation:**

```bash
python validate_projection_scan.py examples/sc_W3_ALocal.xyz \
    examples/sc_W3_ALocal_projection_scan.csv \
    --n-check 500 \
    --seed 2026
```

**Independent 40,000-direction coarse search:**

```bash
python validate_projection_scan.py examples/sc_W3_ALocal.xyz \
    examples/sc_W3_ALocal_projection_scan.csv \
    --n-check 500 \
    --extra-directions 40000 \
    --coarse-points 500
```

The latter is intentionally more computationally expensive.

For the supplied 10,000-direction scan, the 2,000-point recheck agrees in 499
of 500 randomly selected directions. One direction near a projection boundary
changes from `(N,C)=(6,0)` at 1,000 points to `(7,-1)` at 2,000 points. The
independent 40,000-direction, 500-point search still finds minimum `N=3`.

#### Validation options

| Option | Meaning |
|---|---|
| `xyz_file` | Original full-resolution XYZ coordinates. |
| `scan_csv` | CSV produced by `analyze_projection_crossings.py`. |
| `--n-check N` | Randomly select and recompute up to `N` Fibonacci rows using the full-resolution XYZ curve. Default: `500`. |
| `--seed N` | Make the random subset reproducible. Default: `2026`. |
| `--engine` | Select `auto`, `shapely`, or `numpy`, as in the analysis script. |
| `--extra-directions N` | Run an additional independent Fibonacci-sphere search of `N` directions. Default: `0`, which disables it. |
| `--coarse-points N` | Resample size for the optional independent coarse search. Default: `500`. |
| `--progress-every N` | Progress-reporting interval for the optional coarse search. |
| `--output PATH` | Override the default `_projection_scan_validation.txt` report path. |
| `--gui` | Open the graphical interface, pre-filled with any other arguments given. |

A mismatch is not automatically a sign-convention failure. Crossing counts
are piecewise constant in viewing direction and jump at nongeneric boundaries.
A slightly different resampled polyline can move such a boundary enough to
change a count for a direction extremely close to it. Inspect the mismatch
rate, the direction, and whether the global averages/minimum conclusions are
stable.

### `compare_projection_scans.py`

**Function.** Compares two conformations direction by direction. The script
uses only rows whose labels begin with `fib_`, checks that the labels and
viewing normals agree between the two CSV files, and then calculates paired
differences in `N(n)` and `C(n)`. It also draws coordinate-view comparisons,
side-by-side distributions, and Mollweide maps of the direction-matched
differences.

The two scan CSVs should be produced with the same `--n-directions` value.
Because `fibonacci_sphere()` is deterministic, equal scan sizes give the same
`fib_` labels and normals.

```bash
python compare_projection_scans.py \
    examples/sc_W3_ALocal.xyz examples/sc_W3_ALocal_projection_scan.csv \
    examples/sc_W3_R11.xyz examples/sc_W3_R11_projection_scan.csv \
    --label-a "ALocal" \
    --label-b "R11" \
    --analysis-points 1000 \
    --format both \
    --outdir examples/conformation_comparison_figures
```

Output names start with `<label-a>_vs_<label-b>` unless `--stem` is given.
`--outdir` defaults to `conformation_comparison_figures` in the current
directory. The command above writes these files to
`examples/conformation_comparison_figures/`:

- `ALocal_vs_R11_matched_directions.csv` — both curves' counts and
  paired differences for every matched viewing direction.
- `ALocal_vs_R11_comparison_summary.txt` — means, minimum-view
  fractions, exact-match fractions, correlations, and largest differences.
- `ALocal_vs_R11_distribution_comparison.png/.svg` — overlaid signed
  and unsigned distributions.
- `ALocal_vs_R11_direction_difference_maps.png/.svg` — maps of
  `N_R11-N_ALocal` and `C_R11-C_ALocal` on the viewing sphere.
- `ALocal_vs_R11_coordinate_views.png/.svg` — xz and yz projections of
  both conformations.

For the supplied 10,000-direction scans, both conformations have mean signed
crossing sum approximately `+3`, but their unsigned distributions differ
substantially. `ALocal` has mean `N=4.8833` and `N=3` for `14.75%` of sampled
directions; `R11` has mean `N=3.6844` and `N=3` for `62.82%`. Only `38.77%`
of matched directions have the same unsigned count. Thus equal writhe does
not imply equal projection-crossing distributions.

### `writhe_gui.py`

**Function.** Builds the graphical interface for any of the four scripts from
that script's own command-line parser, so the window always offers the same
options as `--help`. Each script calls it when started with no arguments or
with `--gui`; running `python writhe_gui.py` opens all four as tabs. See
[Graphical interface](#graphical-interface).

## Sign convention

The scripts use the Gauss-integral crossing convention

`sign[((t_i x t_j) . (r_i - r_j))]`.

At a projected crossing this is the usual over-strand/under-strand sign in the oriented screen basis, and it is also the sign of the Gauss-writhe integrand. With this convention, `sc_W3_ALocal.xyz` has mean signed crossing sum near `+3`, consistent with direct numerical evaluation of its Gauss writhe. Reversing the arbitrary 2D screen basis must not flip this value. The central result remains: **one projection gives a directional integer `w(n)`, while the 3D writhe is the spherical average of `w(n)` and need not equal any particular projection.**

## Is `sc_W3_ALocal.xyz` a positive supercoil?

Its centerline has positive Gauss writhe (`Wr` approximately `+3`), so it is a **positive-writhe, left-handed plectonemic centerline** under the standard DNA-topology convention. Strictly, an XYZ centerline alone does not contain the ribbon framing or strand linking needed to determine `Tw`, `Lk`, or `Delta Lk`. It can therefore establish positive writhe, while the stronger statement “positively supercoiled DNA” additionally assumes that positive excess linking number is being stored in this writhe.

## Exact CLI sequence for the report graphics

Run these commands from the repository root. They regenerate the files in
`examples/`:

```bash
python analyze_projection_crossings.py examples/sc_W3_ALocal.xyz \
    --n-directions 10000 \
    --analysis-points 1000 \
    --include-axis-views \
    --direction 0.2494,0.8564,0.4521 \
    --direction=-0.2700,0.8643,0.4243 \
    --direction=-0.0244,0.6586,0.7521 \
    --direction=-0.1754,-0.6293,0.7571 \
    --direction 0.1941,0.2748,0.9417

python plot_writhe_projection_views.py examples/sc_W3_ALocal.xyz \
    --scan-csv examples/sc_W3_ALocal_projection_scan.csv \
    --analysis-points 1000 \
    --annotate-crossings \
    --report-panels \
    --outdir examples/sc_W3_ALocal_writhe_figures

python validate_projection_scan.py examples/sc_W3_ALocal.xyz \
    examples/sc_W3_ALocal_projection_scan.csv \
    --n-check 500 \
    --seed 2026 \
    --extra-directions 40000 \
    --coarse-points 500

python analyze_projection_crossings.py examples/sc_W3_R11.xyz \
    --n-directions 10000 \
    --analysis-points 1000 \
    --include-axis-views

python plot_writhe_projection_views.py examples/sc_W3_R11.xyz \
    --scan-csv examples/sc_W3_R11_projection_scan.csv \
    --analysis-points 1000 \
    --annotate-crossings \
    --report-panels \
    --format both \
    --outdir examples/sc_W3_R11_writhe_figures

python validate_projection_scan.py examples/sc_W3_R11.xyz \
    examples/sc_W3_R11_projection_scan.csv \
    --n-check 500 \
    --seed 2026 \
    --extra-directions 40000 \
    --coarse-points 500

python compare_projection_scans.py \
    examples/sc_W3_ALocal.xyz examples/sc_W3_ALocal_projection_scan.csv \
    examples/sc_W3_R11.xyz examples/sc_W3_R11_projection_scan.csv \
    --label-a "ALocal" \
    --label-b "R11" \
    --analysis-points 1000 \
    --format both \
    --outdir examples/conformation_comparison_figures
```

## Negative-writhe example: `sc_L714_Wm2_Logo.xyz`

`examples/sc_L714_Wm2_Logo.xyz` is a shorter closed centerline (2,000 points,
closed length 714 Å) with negative writhe. It shows the same distinctions as
the report examples with the opposite sign. Its outputs were produced with:

```bash
python analyze_projection_crossings.py examples/sc_L714_Wm2_Logo.xyz \
    --n-directions 10000 \
    --analysis-points 1000 \
    --include-axis-views

python plot_writhe_projection_views.py examples/sc_L714_Wm2_Logo.xyz \
    --scan-csv examples/sc_L714_Wm2_Logo_projection_scan.csv \
    --analysis-points 1000 \
    --annotate-crossings \
    --only-custom-views \
    --view xz_reference:0,-1,0 \
    --view yz_reference:1,0,0 \
    --view minimum_crossing:0.8503,-0.0592,0.5229 \
    --view=partial_cancellation:-0.7458,-0.0060,0.6661 \
    --view full_cancellation:0.6340,-0.2019,0.7465 \
    --view opposite_sign:0.2128,-0.2562,0.9429 \
    --format both

python validate_projection_scan.py examples/sc_L714_Wm2_Logo.xyz \
    examples/sc_L714_Wm2_Logo_projection_scan.csv \
    --n-check 500 \
    --seed 2026 \
    --extra-directions 40000 \
    --coarse-points 500
```

The seven built-in views of the plotting script (and therefore
`--report-panels`) were chosen for `sc_W3_ALocal.xyz`, so this example uses
`--only-custom-views` with views picked from its own scan. Every view gives
the same `(N, C)` at 1,000 and 2,000 points. The `xz` view and the four
oblique views also keep their `(N, C)` when the normal is changed by 0.004 in
any one coordinate. The `yz`
view does not: it lies next to a boundary between the `N = 2` and `N = 3`
regions, and a tilt that small can change it to `(3, -3)`.

| View | Normal | `N` | `C` |
|---|---|---|---|
| `xz_reference` | `(0, -1, 0)` | 3 | -3 |
| `yz_reference` | `(1, 0, 0)` | 2 | -2 |
| `minimum_crossing` | `(0.8503, -0.0592, 0.5229)` | 2 | -2 |
| `partial_cancellation` | `(-0.7458, -0.0060, 0.6661)` | 3 | -1 |
| `full_cancellation` | `(0.6340, -0.2019, 0.7465)` | 4 | 0 |
| `opposite_sign` | `(0.2128, -0.2562, 0.9429)` | 5 | +1 |

Results of the 10,000-direction scan:

- mean signed crossing sum `C(n)`: `-2.0003` (the writhe estimate)
- mean unsigned crossing count `N(n)`: `3.8949`
- smallest unsigned count: `2`, seen in `1.19%` of directions
- most common signed sum: `C = -3` (`63.26%` of directions), not `-2`
- `14.51%` of directions give `C = +1`, the opposite sign to `Wr`
- looking along the superhelix axis (`xy` view, `+z` normal) gives
  `N = 177`, `C = +1`

The full-resolution recheck agrees in all 500 checked directions, and the
independent 40,000-direction search also finds minimum `N = 2`.

## Recommended workflow

1. Run `analyze_projection_crossings.py` at 1000 points for the main spherical scan.
2. Run `plot_writhe_projection_views.py` to generate the teaching figures.
3. Run `validate_projection_scan.py` against the original 2000-point curve to confirm the reduced-resolution scan.
4. For two conformations, generate both scans with the same number of
   directions and run `compare_projection_scans.py`.
5. If needed for a manuscript/reviewer, add the independent large coarse search.

## Tests

Run all tests from the repository root:

```bash
python -m unittest -v
```

`test_crossing_signs.py` checks the corrected positive signs in the xz and yz
projections, invariance under reversing the viewing normal, the six-crossing
cancellation view, and a positive spherical mean. `test_gui.py` checks that
the graphical interface offers every command-line option, builds commands
that parse back to the intended values, and pre-fills the form from `--gui`
arguments. When a display is available, it also runs and stops a short
analysis from a hidden window.

## Performance and accuracy guidance

- Runtime grows strongly with the number of curve segments and approximately
  linearly with the number of viewing directions.
- Begin with `--n-directions 1000 --analysis-points 500` when testing a new
  file or command.
- Use `--n-directions 10000 --analysis-points 1000` for the supplied report
  graphics and statistics.
- Use `validate_projection_scan.py` to compare a reduced scan against the
  original coordinates instead of assuming resampling has no effect.
- Shapely is normally much faster for large jobs. The NumPy fallback is useful
  for portability and testing.
- Do not interpret the unsigned mean as writhe. The mean of `C(n)` estimates
  `Wr`; the mean of `N(n)` estimates average crossing number.

## Common command-line issues

- **A negative vector is treated as an option:** write
  `--direction=-0.27,0.86,0.42` instead of putting a space after
  `--direction`.
- **No histograms or maps appear:** supply `--scan-csv` to the plotting
  script and confirm the CSV path is correct.
- **Standard views are regenerated:** add `--only-custom-views`. This flag is
  optional and affects projection images only.
- **`--report-panels` fails with missing built-in views:** remove
  `--only-custom-views`; the combined report panels require the seven standard
  normals.
- **Changing plotting `--analysis-points` does not change a histogram:** this
  is expected. Histogram/map values come from the CSV. Regenerate the CSV with
  the analysis script to change those values.
- **A custom output filename is needed:** use `--csv-out`, `--summary-out`,
  `--output`, or `--outdir`, depending on the script.
- **A window opened instead of the script running:** the script was started
  without arguments or with `--gui`. Give the arguments on the command line
  (and leave out `--gui`) to run it directly.
- **The window does not open:** the script prints the reason followed by its
  `--help` text. Install `tkinter` (see [Installation](#installation)) or use
  the command line.

## License

Released under the MIT License; see [`LICENSE`](LICENSE).
