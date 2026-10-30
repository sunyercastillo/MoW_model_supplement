# Mo–W systematics of subduction-zone fluids — forward models

Supplementary material to *<Article title>* (Geochimica et Cosmochimica Acta,
2026). The repository contains the two forward models used in the article, the
Excel workbook that reproduces every number of the sensitivity grids, the
DataGraph tables behind the figure, and an interactive dashboard to explore the
models with other inputs.

**Live dashboard (no installation):** `<Streamlit Community Cloud URL>`

## The two models

1. **Mo-W Thermodynamic model.** Fluid–mineral partitioning of Bali et al.
   (2011, 2012) at 1000 °C and 26.1 kbar for an eclogite residue
   (clinopyroxene + garnet + rutile). Three sensitivity grids in Mo/W space, each
   drawn as ΔFMQ bands (−2 to +3) crossed by constant-parameter isopleths:
   residual rutile (x = Nb/La), fluid salinity (x = U/Th) and sediment addition
   to the N-MORB protolith (x = U/Th) with a selectable carbonate:hemipelagic
   proportion.
2. **Two-stage fractionation model.** Rayleigh dehydration of a slab segment
   (Stage 1, fractionation factor α) followed by re-equilibration of the fluid in a
   sediment-buffered redox filter (Stage 2, α_r), tracking Mo and δ⁹⁸/⁹⁵Mo of the
   cumulative fluid.

Both models draw the same lava data (Nicaragua, Costa Rica) from `data/arc_data.xlsx`.

## Repository layout

```
dashboard.py                  Streamlit application (entry point)
launch_dashboard.sh / .bat    one-click local launchers (macOS/Linux, Windows)
requirements.txt              Python dependencies
thermodynamic_model/
    model.py                  equations, constants and grid builders of the Mo-W Thermodynamic model
    build_workbook.py         writes outputs/ (workbook, DataGraph tables, validation report)
two_stage_model/
    bridge.py                 reads data/ and feeds the two-stage model
    fluid_fractionation.py    two-stage model and its plot
    CORE_scientific_calculations_engine.py, sediment_endmembers.py, arc_data_overlay.py
data/                         editable inputs (see "Changing inputs")
    model_inputs.xlsx         constants, end-members, grid axes
    sample_input_data.xlsx    slab segments for the two-stage model
    arc_data.xlsx             lava data and plotting styles
outputs/
    MoW_thermodynamic_model_grids.xlsx   workbook with live formulas (Inputs, Rutile, Salinity,
                                         Pure_carbonate, C60_H40, DG_*, Checks)
    datagraph/DG_*.txt                   tab-separated tables for DataGraph, one per grid
    validation_report.md                 checks performed on grids and workbook
docs/
    audit_inputs_and_equations.md        equations, constants, and inconsistencies found in the sources
    DataGraph_instructions.md            how the figure is assembled from the tables
.streamlit/config.toml        dashboard theme
CITATION.cff · LICENSE (CC BY 4.0)
```

`thermodynamic_model/model.py` is the single source of the Mo-W Thermodynamic
model: the dashboard, the workbook builder and the DataGraph tables all use it,
and the Excel formulas were checked against it (see the validation report).

## Running the dashboard locally

Requirements: Python 3.9 or newer (https://www.python.org/downloads/). Nothing
else; the dashboard runs in the web browser.

- macOS / Linux: `./launch_dashboard.sh`
- Windows: double-click `launch_dashboard.bat`

The launcher creates a local virtual environment (`.venv`) on first use,
installs `requirements.txt` (one-off, needs internet) and opens
http://localhost:8501. Unpack the folder in a path without spaces if you want
the PNG/SVG export buttons of the two-stage model to work (the camera icon of
every plot works regardless). By hand:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/streamlit run dashboard.py
```

## Reproducing the published grids

```bash
.venv/bin/python -m thermodynamic_model.build_workbook
```

rebuilds `outputs/` from `data/model_inputs.xlsx`: the workbook with live,
cell-referenced formulas (every constant is a named cell on `Inputs`), the four
DataGraph tables and `outputs/validation_report.md`. The report lists the
structural and geometric checks on the grids, an independent evaluation of every
workbook formula against the Python model, and (when Microsoft Excel is
available) a recalculation of the workbook by Excel itself. The dashboard's
*Plotted XY pairs (TSV)* downloads are identical to the files in
`outputs/datagraph/`.

## Changing inputs

All inputs are read from the three workbooks in `data/`. Save the workbook and
reload the browser page; no Python editing is needed.

### Model constants and end-members — `data/model_inputs.xlsx`

Sheet `Parameters` (`Name`, `Value`, `Unit`, `Note`); change column B.

| Want to … | Edit |
|---|---|
| Use 1.5 % or 2 % slab fluid in the source instead of 1 % | `f_fluid` → 0.015 or 0.02 |
| Change the reference salinity or rutile content | `NaCl_wt`, `Rut_wt` |
| Change the cpx:grt ratio of the silicate assemblage | `cpx_share` (0.70 = 70:30) |
| Use the Appendix C Mo coefficient | `Mo_a` → 0.435 |
| Use a batch fluid fraction other than 0 | `F_batch` |
| Change temperature or pressure | `T_C`, `P_bar` |

Sheet `Endmembers` (µg/g): N-MORB, depleted mantle, carbonate and hemipelagic
sediment (Mo, W, U, Th, Nb, La). Replace a column to model another sediment or
source at the same conditions; the dashboard slider mixes the two sediment
columns.

Sheet `Grids`: ΔFMQ values, sediment additions, rutile contents and salinities
of the grids. Values are read left to right; blank cells shorten a list; the
sediment list must contain 0.

The build script never overwrites `data/model_inputs.xlsx`; delete it to
restore the published defaults (it is recreated on the next run). The `Checks`
sheet and the ΔFMQ 0 comparison in the report assume the published grid axes.
The workbook in `outputs/` can also be edited on its own through its `Inputs`
sheet.

### Lava data — `data/arc_data.xlsx`

Sheet `ArcData`, one row per sample: `SampleID`, `Arc`, `Group`, `Mo_ppm`,
`W_ppm`, `U_ppm`, `Th_ppm`, `Nb_ppm`, `La_ppm`, and for the two-stage model
`d98Mo`, `d98Mo_2SE`. Rows without Mo or W are skipped. A new value in `Arc`
becomes a new legend entry; Nicaragua and Costa Rica use the colours of the
article, other arcs take a `Style` row with `KeyType = Arc`, `Key = <Arc>`,
`Color`, `OutlineColor`, or a default palette. The `Opacity` of the `Style`
rows with `KeyType = Group` sets the transparency of that group's samples in
all plots (1 = opaque, 0.5 = half-transparent, as in the article).

### Slab segments — `data/sample_input_data.xlsx`

Sheet `Part2_Segments`, one row per slab segment: `segment`, `P_bar`, `T_C`
(`T_K` and `phi` are recomputed if Excel has not cached the formulas),
`ΔFMQ (log units)`, `ΔFMQ_filter`, `C0_MO_filter`, `mode_grt`, `mode_cpx`,
`mode_rut`, `C0_MO` … `C0_La`, `d98Mo`, `d98Mo 2SD`, `rock_id`. A row is used
when pressure, temperature, ΔFMQ, the three modes and C0_MO, C0_Ce are filled;
the segment is selected with the tick box in the dashboard. Sheet `Part1`: the
first complete row provides the reference partition-coefficient table used as
fallback by the two-stage model.

## Fixed choices worth knowing

- Mo fO₂ coefficient 0.44 (article), not the 0.435 of Bali et al. (2012)
  Appendix C; `docs/audit_inputs_and_equations.md` §2.1 quantifies the effect
  and `Inputs!Mo_a` lets the reader test either value.
- The residual-rutile panel is drawn with the slab-fluid composition; the
  salinity and sediment panels with the metasomatised source (1 % fluid + 99 %
  depleted mantle).
- End-member concentrations are the Table B1 values of the article.

## Citation and licence

Please cite the article and this repository (`CITATION.cff`). Content is
released under the Creative Commons Attribution 4.0 International licence
(`LICENSE`).

## References

Bali E., Audétat A., Keppler H. (2011) The mobility of U and Th in subduction
zone fluids: an indicator of oxygen fugacity and fluid salinity. *Contrib.
Mineral. Petrol.* 161, 597–613.

Bali E., Keppler H., Audétat A. (2012) Molybdenum and tungsten solubility in
aqueous fluids at high pressure and temperature and implications for arc
magmatism. *Earth Planet. Sci. Lett.* 351–352, 195–207.
