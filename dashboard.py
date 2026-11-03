#!/usr/bin/env python3
"""
Interactive supplement: Mo–W systematics of subduction-zone fluids.

Mo-W Thermodynamic model (Bali et al. 2011, 2012 fluid–mineral partitioning):
    three sensitivity grids in Mo/W space — residual rutile, fluid salinity and
    sediment addition — drawn as ΔFMQ bands with constant-parameter isopleths.
Two-stage fractionation model: Rayleigh slab dehydration followed by a redox filter.

Run:  streamlit run dashboard.py   (or ./launch_dashboard.sh)
"""

import os
import math

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from thermodynamic_model.model import (Model, MODEL_INPUTS_XLSX, param_tag, dg_table, tsv_text, grid_dataframe_rows,
                                       write_model_inputs_template)

HERE = os.path.dirname(os.path.abspath(__file__))
ARC_XLSX = os.path.join(HERE, "data", "lava_data.xlsx")

# Colours sampled from the published figure
BAND_COLORS = ["#e3edf6", "#eef5fa", "#fafafa", "#fcede8", "#f9dbd5"]   # ΔFMQ −2→−1 … +2→+3
CURVE_COLOR = "#4d4d4d"
ISO_COLOR = "#808080"
REF_COLOR = "#9678b8"
ARC_STYLE = {   # colours of the published figure; other arcs take the Style sheet (KeyType = Arc) or the palette
    "Nicaragua": {"fill": "#fca9ff", "line": "#9c329a"},
    "CostaRica": {"fill": "#93e65d", "line": "#3a7a1e"},
}
ARC_LABEL = {"Nicaragua": "Nicaragua", "CostaRica": "Costa Rica"}
ARC_PALETTE = ["#f0c050", "#60c0e0", "#c0a0e0", "#e08060", "#a0a0a0"]
PLOT_PX = 480   # fixed panel size in pixels (does not follow the window width)

st.set_page_config(page_title="Mo–W systematics — forward models", layout="wide")


# ──────────────────────────────────────────────────────────────────────────────
# Data
# ──────────────────────────────────────────────────────────────────────────────
@st.cache_data
def load_arc_data(_mtime):
    """ArcData sheet: one row per lava sample (SampleID, Arc, Mo_ppm, W_ppm, U_ppm, Th_ppm, Nb_ppm, La_ppm).
    Colours: ARC_STYLE for the two arcs of the paper, then rows of the Style sheet with KeyType = 'Arc',
    then a default palette."""
    df = pd.read_excel(ARC_XLSX, sheet_name="ArcData")
    for c in ("Mo_ppm", "W_ppm", "U_ppm", "Th_ppm", "Nb_ppm", "La_ppm"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["MoW"] = df["Mo_ppm"] / df["W_ppm"]
    df["UTh"] = df["U_ppm"] / df["Th_ppm"]
    df["NbLa"] = df["Nb_ppm"] / df["La_ppm"]
    df = df.dropna(subset=["MoW", "Arc"])
    df["opacity"] = 1.0
    styles = dict(ARC_STYLE)
    try:
        sty = pd.read_excel(ARC_XLSX, sheet_name="Style")
        for _, r in sty[sty["KeyType"].astype(str).str.lower() == "arc"].iterrows():
            styles.setdefault(str(r["Key"]), {"fill": str(r["Color"]), "line": str(r.get("OutlineColor", "#333333"))})
        # per-Group opacity (e.g. 0.5 for samples drawn half-transparent in the paper)
        grp = sty[sty["KeyType"].astype(str).str.lower() == "group"].set_index("Key")["Opacity"].astype(float)
        if "Group" in df.columns:
            df["opacity"] = df["Group"].map(grp).fillna(1.0)
    except Exception:
        pass
    k = 0
    for arc in df["Arc"].unique():
        if arc not in styles:
            styles[arc] = {"fill": ARC_PALETTE[k % len(ARC_PALETTE)], "line": "#555555"}
            k += 1
    return df, styles


@st.cache_resource
def model(_mtime):
    return Model()


def dg_headers(param_name, xname, yname, iso_prefix, iso_vals):
    h = ["FMQ", param_name, xname, yname, "Band_mask"]
    for b in range(1, 6):
        h += [f"Band{b}_X", f"Band{b}_Y"]
    for v in iso_vals:
        h += [f"{iso_prefix}{param_tag(v)}_X", f"{iso_prefix}{param_tag(v)}_Y"]
    return h


# ──────────────────────────────────────────────────────────────────────────────
# Figure
# ──────────────────────────────────────────────────────────────────────────────
def grid_figure(cols, param_vals, dfmq_values, iso_vals, ref_value, iso_label, x_title, y_title, x_range, y_range,
                panel_label, caption, arc_df=None, arc_styles=None, arc_x="UTh", show_legend=False):
    """Bands → solid ΔFMQ curves → dashed isopleths → data, as in the published figure."""
    d = dict(cols)
    n = len(param_vals)
    fig = go.Figure()
    for b in range(1, len(dfmq_values)):
        fig.add_trace(go.Scatter(x=d[f"Band{b}_X"], y=d[f"Band{b}_Y"], mode="lines", fill="toself",
                                 fillcolor=BAND_COLORS[(b - 1) % len(BAND_COLORS)], line=dict(width=0), hoverinfo="skip",
                                 showlegend=False))
    for i, dfmq in enumerate(dfmq_values):
        xs, ys = d["X"][i * n:(i + 1) * n], d["MoW"][i * n:(i + 1) * n]
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=CURVE_COLOR, width=1.3),
                                 hovertemplate=f"ΔFMQ {dfmq:+g}<br>%{{x:.3g}}, %{{y:.3g}}<extra></extra>",
                                 showlegend=False))
        # ΔFMQ label just outside the grid, on the side where the curve ends
        left_end = xs[-1] < xs[0]
        fig.add_annotation(x=math.log10(xs[-1]), y=math.log10(ys[-1]), text=f"{dfmq:+g}" if dfmq else "0",
                           showarrow=False, xanchor="right" if left_end else "left", xshift=-6 if left_end else 6,
                           font=dict(size=11, color=CURVE_COLOR))
    for k, v in enumerate(iso_vals):
        xs, ys = d[f"ISO{param_tag(v)}_X"], d[f"ISO{param_tag(v)}_Y"]
        is_ref = v == ref_value
        col = REF_COLOR if is_ref else ISO_COLOR
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=col, width=1.4 if is_ref else 0.9, dash="dash"),
                                 hovertemplate=f"{iso_label(v)}<br>%{{x:.3g}}, %{{y:.3g}}<extra></extra>", showlegend=False))
        fig.add_trace(go.Scatter(x=[xs[-1]], y=[ys[-1]], mode="markers", hoverinfo="skip", showlegend=False,
                                 marker=dict(size=6, color=col, line=dict(color="white", width=0.5))))
        fig.add_annotation(x=math.log10(xs[-1]), y=math.log10(ys[-1]), text=iso_label(v), showarrow=False,
                           yanchor="bottom", yshift=4 + 12 * (k % 2),  # staggered to avoid overlap
                           font=dict(size=10, color=REF_COLOR if is_ref else ISO_COLOR))
    if arc_df is not None:
        for arc in arc_df["Arc"].unique():
            sty = arc_styles[arc]
            sub = arc_df[arc_df["Arc"] == arc].dropna(subset=[arc_x])
            fig.add_trace(go.Scatter(x=sub[arc_x], y=sub["MoW"], mode="markers", name=ARC_LABEL.get(arc, arc),
                                     marker=dict(size=9, color=sty["fill"], opacity=sub["opacity"].tolist(),
                                                 line=dict(color=sty["line"], width=1)),
                                     text=sub["SampleID"], hovertemplate="%{text}<br>%{x:.3g}, %{y:.3g}<extra></extra>"))
    fig.add_annotation(xref="paper", yref="paper", x=0.02, y=0.98, text=f"<b>{panel_label}</b>", showarrow=False,
                       font=dict(size=14), bordercolor="#333", borderwidth=1, borderpad=3, bgcolor="white")
    fig.add_annotation(xref="paper", yref="paper", x=0.98, y=0.03, text=f"<b>{caption}</b>", showarrow=False,
                       xanchor="right", font=dict(size=12))
    fig.update_layout(
        template="simple_white", width=PLOT_PX, height=PLOT_PX, margin=dict(l=60, r=30, t=30, b=55),
        paper_bgcolor="white", plot_bgcolor="white", font=dict(color="#222222"),
        showlegend=show_legend,
        legend=dict(x=0.03, y=0.90, xanchor="left", yanchor="top", bgcolor="rgba(255,255,255,0.8)", font=dict(size=10)),
        xaxis=dict(type="log", title=x_title, range=[math.log10(x_range[0]), math.log10(x_range[1])], dtick=1,
                   tickformat=".3~g", showgrid=False, showline=True, mirror=True, linecolor="black", linewidth=1,
                   ticks="outside", minor=dict(ticks="outside", showgrid=False)),
        yaxis=dict(type="log", title=y_title, range=[math.log10(y_range[0]), math.log10(y_range[1])], dtick=1,
                   tickformat=".3~g", showgrid=False, showline=True, mirror=True, linecolor="black", linewidth=1,
                   ticks="outside", minor=dict(ticks="outside", showgrid=False)),
    )
    return fig


def panel(column, fig, headers, cols, grid_rows, param_label, file_stub, key):
    with column:
        st.plotly_chart(fig, use_container_width=False, theme=None,
                        config={"toImageButtonOptions": {"format": "svg", "scale": 2}, "displaylogo": False})
        c1, c2 = st.columns(2)
        c1.download_button("Plotted XY pairs (TSV)", tsv_text(headers, cols), file_name=f"{file_stub}.txt",
                           mime="text/tab-separated-values", key=f"dl_xy_{key}",
                           help="FMQ curves, Band1–5 polygons with Band_mask, and the dashed isopleths — "
                                "the DataGraph layout used for the article figure.")
        full = pd.DataFrame(grid_dataframe_rows(grid_rows, param_label))
        c2.download_button("Full calculation table (CSV)", full.to_csv(index=False), file_name=f"{file_stub}_full.csv",
                           mime="text/csv", key=f"dl_full_{key}")


# ──────────────────────────────────────────────────────────────────────────────
# Page
# ──────────────────────────────────────────────────────────────────────────────
st.title("Mo-W systematics of subduction-zone fluids")

write_model_inputs_template()   # restores data/model_inputs.xlsx from the defaults if it was deleted
m = model(os.path.getmtime(MODEL_INPUTS_XLSX))
c = m.c
arc_mtime = os.path.getmtime(ARC_XLSX) if os.path.exists(ARC_XLSX) else 0

with st.expander("Mo-W Thermodynamic model", expanded=True):
    st.markdown(
        f"Fluid–mineral partitioning of Bali et al. (2011, 2012). Fixed conditions: T = {c['T_C']:g} °C, "
        f"P = {c['P_bar']:g} bar, cpx:grt = {c['cpx_share']:.0%}:{1 - c['cpx_share']:.0%} of the silicate assemblage, "
        f"{c['f_fluid']:.0%} slab fluid + {1 - c['f_fluid']:.0%} depleted mantle in the source. "
        f"Reference case: {c['NaCl_wt']:g} wt.% NaCl, {c['Rut_wt']:g} wt.% rutile, N-MORB protolith "
        "(purple dashed line in each panel). Bands span ΔFMQ −2 to +3 in one-log-unit steps; "
        "dashed lines hold the varied parameter constant. Panel A plots the slab-fluid composition itself; "
        "panels B and C plot the metasomatised source."
    )
    ctl1, ctl2, ctl3 = st.columns([1, 2, 1])
    show_data = ctl1.checkbox("Show lava data (Nicaragua, Costa Rica)", value=True)
    arc_df, arc_styles = load_arc_data(arc_mtime) if (show_data and os.path.exists(ARC_XLSX)) else (None, None)
    preset = ctl2.radio("Sediment end-member (panel C)",
                        ["Pure carbonate (100C:0H)", "60C:40H", "Pure hemipelagic (0C:100H)", "Custom"],
                        horizontal=True, index=0, key="sed_preset")
    if preset.startswith("Pure carbonate"):
        carb_pct = 100
    elif preset == "60C:40H":
        carb_pct = 60
    elif preset.startswith("Pure hemipelagic"):
        carb_pct = 0
    else:
        carb_pct = ctl3.slider("Carbonate proportion of the sediment (%)", 0, 100, 60, 5, key="sed_carb_pct")

    colA, colB, colC = st.columns(3)

    # A — residual rutile (fluid composition)
    grid_rut = m.rutile_grid()
    cols_rut = dg_table(grid_rut, "Fluid_NbLa", m.rut_wt, m.rut_wt, ykey="Fluid_MoW")
    hdr_rut = dg_headers("Rutile_wt_pct", "Fluid_NbLa", "Fluid_MoW", "Rut", m.rut_wt)
    fig_a = grid_figure(cols_rut, m.rut_wt, m.dfmq, m.rut_wt, c["Rut_wt"],
                        lambda v: f"{v:g}%" if v != c["Rut_wt"] else f"{v:g}% Rut",
                        "Nb/La (fluid)", "Mo/W (fluid)", (0.07, 7), (0.5, 300), "A", "Residual rutile",
                        arc_df, arc_styles, arc_x="NbLa")
    panel(colA, fig_a, hdr_rut, cols_rut, grid_rut, "Rutile (wt.%)", "A_rutile_grid", "A")

    # B — fluid salinity
    grid_sal = m.salinity_grid()
    cols_sal = dg_table(grid_sal, "UTh", m.nacl_wt, m.nacl_wt)
    hdr_sal = dg_headers("NaCl_wt_pct", "UTh", "MoW", "NaCl", m.nacl_wt)
    fig_b = grid_figure(cols_sal, m.nacl_wt, m.dfmq, m.nacl_wt, c["NaCl_wt"],
                        lambda v: f"{v:g}%" if v != c["NaCl_wt"] else f"NaCl {v:g}%",
                        "U/Th", "Mo/W", (0.3, 50), (0.5, 300), "B", "Fluid salinity", arc_df, arc_styles)
    panel(colB, fig_b, hdr_sal, cols_sal, grid_sal, "NaCl (wt.%)", "B_salinity_grid", "B")

    # C — sediment addition, carbonate:hemipelagic proportion selectable
    sediment = m.mixed_sediment(carb_pct / 100.0)
    grid_sed = m.sediment_grid(sediment)
    sed_iso = [v for v in m.sed_wt if v != 0]
    cols_sed = dg_table(grid_sed, "UTh", sed_iso, m.sed_wt)
    hdr_sed = dg_headers("Sediment_wt_pct", "UTh", "MoW", "Sed", sed_iso)
    caption = f"Sediment addition ({carb_pct:g}C:{100 - carb_pct:g}H)"
    fig_c = grid_figure(cols_sed, m.sed_wt, m.dfmq, sed_iso, None, lambda v: f"{v:g}", "U/Th", "Mo/W", (0.3, 50), (0.5, 300),
                        "C", caption, arc_df, arc_styles, show_legend=True)
    # 0 wt.% reference line (first point of every ΔFMQ curve) drawn in purple, as in the figure
    x0 = [r["UTh"] for d_, w, r in grid_sed if w == 0]
    y0 = [r["MoW"] for d_, w, r in grid_sed if w == 0]
    fig_c.add_trace(go.Scatter(x=x0, y=y0, mode="lines", line=dict(color=REF_COLOR, width=1.4, dash="dash"),
                               hovertemplate="0 wt.% sediment<br>%{x:.3g}, %{y:.3g}<extra></extra>", showlegend=False))
    fig_c.add_trace(go.Scatter(x=[x0[-1]], y=[y0[-1]], mode="markers", hoverinfo="skip", showlegend=False,
                               marker=dict(size=6, color=REF_COLOR, line=dict(color="white", width=0.5))))
    fig_c.add_annotation(x=math.log10(x0[-1]), y=math.log10(y0[-1]), text="0", showarrow=False, yanchor="bottom",
                         yshift=4, font=dict(size=10, color=REF_COLOR))
    panel(colC, fig_c, hdr_sed, cols_sed, grid_sed, "Sediment (wt.%)", f"C_sediment_grid_{carb_pct:g}C{100 - carb_pct:g}H", "C")

    with st.expander("End-member compositions and sediment mixture used (µg/g)"):
        em = pd.DataFrame(m.em).rename(columns={"NMORB": "N-MORB", "DM": "Depleted mantle", "Carb": "Carbonate",
                                                "Hemi": "Hemipelagic"})
        em[f"Sediment {carb_pct:g}C:{100 - carb_pct:g}H"] = pd.Series(sediment)
        st.dataframe(em.loc[["Mo", "W", "U", "Th", "Nb", "La"]].style.format("{:.4g}"), use_container_width=True)
        st.caption("Nb and La enter only panel A (pure N-MORB protolith); the sediment Nb and La values are not used.")
        if m.inputs_file:
            st.caption(f"Constants, end-members and grid axes are read from `{os.path.relpath(m.inputs_file, HERE)}` — "
                       "edit that workbook (e.g. f_fluid, sediment compositions) and reload the page.")

with st.expander("Two-stage fractionation model", expanded=False):
    try:
        from two_stage_model.bridge import render_two_stage_model
        render_two_stage_model()
    except Exception as e:  # keep the first model usable if the second fails
        st.error(f"Two-stage fractionation model could not be rendered: {e}")
