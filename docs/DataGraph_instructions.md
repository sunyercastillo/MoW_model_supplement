# DataGraph construction of the Mo/W sensitivity grids

Four tab-separated tables are provided in `outputs/datagraph/`, one per grid. Each is
copied as a block into a DataGraph data table (File → Import, or paste into an
empty table with the first row as column names). Column lengths differ by
design: the main curves occupy 30–42 rows, the band polygons 11–15 rows and the
isopleths 6 rows; the remaining cells are empty and must stay empty.

| File | Grid | X column | Y column | Parameter column | Isopleth prefix |
|---|---|---|---|---|---|
| `DG_Rutile.txt` | Residual rutile (panel A) | `Fluid_NbLa` | `Fluid_MoW` | `Rutile_wt_pct` | `Rut` (0, 0_25, 0_5, 1, 2_5, 3_5, 5) |
| `DG_Salinity.txt` | Fluid salinity (panel B) | `UTh` | `MoW` | `NaCl_wt_pct` | `NaCl` (1, 3, 5, 10, 15, 20) |
| `DG_Pure_carbonate.txt` | Sediment, pure carbonate (panel C) | `UTh` | `MoW` | `Sediment_wt_pct` | `Sed` (1, 2_5, 5, 10) |
| `DG_C60_H40.txt` | Sediment, 60 % carbonate + 40 % hemipelagic (panel D) | `UTh` | `MoW` | `Sediment_wt_pct` | `Sed` (1, 2_5, 5, 10) |

Panel A is plotted with the slab-fluid composition (no mantle mixing), as in
the published figure; panels B–D are plotted with the metasomatised source
(1 % fluid + 99 % depleted mantle). The column names carry this distinction.

All axes are logarithmic.

## 1. Coloured ΔFMQ bands (draw first, so they sit behind the lines)

Five *Plot* commands, one per band, using `Band1_X`/`Band1_Y` … `Band5_X`/`Band5_Y`.

- Mask: `Band_mask` equals 1 (the `Band_mask` column holds 1 for the polygon rows and nothing afterwards). Do not mask with the `FMQ` or the parameter column.
- Line style: connect every point, in table order. Do **not** enable "connect only while X increases" and do **not** sort by X — the row order is the polygon (lower ΔFMQ boundary forward through the parameter, upper boundary backwards, first point repeated to close).
- Fill: "fill inside the path" (closed polygon). No symbols. No outline, or an outline of 0.25 pt in the fill colour.
- Fill colours sampled from the published figure:

| Band | Between ΔFMQ | Fill |
|---|---|---|
| Band1 | −2 and −1 | `#e3edf6` |
| Band2 | −1 and 0 | `#eef5fa` |
| Band3 | 0 and +1 | `#fafafa` |
| Band4 | +1 and +2 | `#fcede8` |
| Band5 | +2 and +3 | `#f9dbd5` |

Each sediment band has 11 points (5 + 5 + 1); the salinity bands have 13 and the rutile bands 15, following the same construction.

## 2. Solid ΔFMQ curves

One *Plot* command with X = `UTh` (or `Fluid_NbLa`), Y = `MoW` (or `Fluid_MoW`).

- Split by the `FMQ` column ("Group by" / "Split by" → `FMQ`) so that each ΔFMQ value becomes its own curve; or use six masked plots `FMQ = −2`, …, `FMQ = 3`.
- Straight solid line, 0.75 pt, colour `#4d4d4d`. Connect only within a ΔFMQ series. The rows of each series are already ordered by increasing parameter (0 → 1 → 2.5 → 5 → 10 wt.% sediment; 0 → 5 wt.% rutile; 1 → 20 wt.% NaCl).
- Label the right-hand end of each curve with its ΔFMQ value (−2, −1, 0, +1, +2, +3).

## 3. Dashed isopleths (constant parameter, ΔFMQ −2 → +3)

One *Plot* command per isopleth pair, e.g. `Sed2_5_X` / `Sed2_5_Y`.

- Connect every point (six rows, ordered ΔFMQ −2 → +3). No fill, no symbols.
- Dashed grey line, 0.5 pt, colour `#808080`.
- The reference case of each grid is drawn as a dashed **purple** line (`#9678b8`): `Rut2_5` in panel A, `NaCl10` in panel B. In panels C and D the reference case is 0 wt.% sediment; no `Sed0` columns are provided — the 0 wt.% points are the first point of every ΔFMQ curve and the left edge of Band1–Band5. To draw the purple reference line in C/D, plot the main `UTh`/`MoW` columns with the mask `Sediment_wt_pct = 0`.
- Label the upper end of each isopleth with its parameter value (e.g. "2.5 % Rut", "10 %", "2.5").

## 4. Data overlay

Lava data (Nicaragua and Costa Rica) are plotted last as filled circles, 6 pt:
Nicaragua fill `#fca9ff`, outline `#9c329a`; Costa Rica fill `#93e65d`, outline
`#3a7a1e`; outline 0.5 pt.

## 5. Axis ranges used in the published figure

| Panel | X | Y |
|---|---|---|
| A | Nb/La 0.07 – 7 | Mo/W 0.5 – 300 |
| B | U/Th 0.3 – 50 | Mo/W 0.5 – 300 |
| C, D | U/Th 0.3 – 50 | Mo/W 0.5 – 300 |

## 6. Regenerating the tables

`python -m thermodynamic_model.build_workbook` (run from the repository root)
rebuilds the workbook, the four tables and `outputs/validation_report.md` from
`data/model_inputs.xlsx` (defaults in `thermodynamic_model/model.py`). Every
value in the tables equals the corresponding cell of the `DG_*` sheets of
`outputs/MoW_thermodynamic_model_grids.xlsx`.
