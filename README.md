# Mo–W systematics of subduction-zone fluids — forward models

Supplementary material to *<Article title>* (Geochimica et Cosmochimica Acta,
2026): the two forward models used in the article as an interactive dashboard
that can be run with other inputs. The Excel workbook reproducing the
sensitivity grids cell by cell, and the DataGraph tables of the figure, are
supplied as a separate supplementary file.

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

Both models draw the same lava data (Nicaragua, Costa Rica) from `data/lava_data.xlsx`.

## Repository layout

```
dashboard.py                  Streamlit application (entry point)
launch_dashboard.sh / .bat    one-click local launchers (macOS/Linux, Windows)
requirements.txt              Python dependencies
thermodynamic_model/
    model.py                  equations, constants and grid builders of the Mo-W Thermodynamic model
two_stage_model/
    bridge.py                 reads data/ and renders the two-stage model
    fractionation.py          Rayleigh dehydration, redox filter and the δ⁹⁸/⁹⁵Mo plot
    partitioning_engine.py    Bali et al. (2011, 2012) partition coefficients for each slab segment
    sediment_endmembers.py    sediment end-members and mixing lines
    lava_data.py              lava data overlay
data/                         editable inputs (see "Changing inputs")
    model_inputs.xlsx         constants, end-members, grid axes
    slab_segments.xlsx    slab segments for the two-stage model
    lava_data.xlsx             lava data and plotting styles
docs/audit_inputs_and_equations.md   equations, constants, and inconsistencies found in the sources
.streamlit/config.toml        dashboard theme
CITATION.cff · LICENSE (CC BY 4.0)
```

`thermodynamic_model/model.py` holds the Mo-W Thermodynamic model; the Excel
workbook of the separate supplementary file was generated from, and verified
against, the same definitions.

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

Delete `data/model_inputs.xlsx` to restore the published defaults; the
dashboard recreates it on the next start.

### Lava data — `data/lava_data.xlsx`

Sheet `ArcData`, one row per sample: `SampleID`, `Arc`, `Group`, `Mo_ppm`,
`W_ppm`, `U_ppm`, `Th_ppm`, `Nb_ppm`, `La_ppm`, and for the two-stage model
`d98Mo`, `d98Mo_2SE`. Rows without Mo or W are skipped. A new value in `Arc`
becomes a new legend entry; Nicaragua and Costa Rica use the colours of the
article, other arcs take a `Style` row with `KeyType = Arc`, `Key = <Arc>`,
`Color`, `OutlineColor`, or a default palette. The `Opacity` of the `Style`
rows with `KeyType = Group` sets the transparency of that group's samples in
all plots (1 = opaque, 0.5 = half-transparent, as in the article).

### Slab segments — `data/slab_segments.xlsx`

Sheet `Part2_Segments`, one row per slab segment: `segment`, `P_bar`, `T_C`
(`T_K` and `phi` are recomputed if Excel has not cached the formulas),
`ΔFMQ (log units)`, `ΔFMQ_filter`, `C0_MO_filter`, `mode_grt`, `mode_cpx`,
`mode_rut`, `C0_MO` … `C0_La`, `d98Mo`, `d98Mo 2SD`, `rock_id`. A row is used
when pressure, temperature, ΔFMQ, the three modes and C0_MO, C0_Ce are filled;
the segment is selected with the tick box in the dashboard. Sheet `Part1`: the
first complete row provides the reference partition-coefficient table used as
fallback by the two-stage model.

## Downloads from the dashboard

Under each grid: *Plotted XY pairs (TSV)* gives the exact plotted coordinates
in the DataGraph layout (ΔFMQ curves, Band1–5 polygons with `Band_mask`, and
the dashed isopleths); *Full calculation table (CSV)* lists every intermediate
quantity of every grid point. Plotly's camera icon exports each panel as SVG.
The two-stage model has its own copy block for the plotted lines.

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
