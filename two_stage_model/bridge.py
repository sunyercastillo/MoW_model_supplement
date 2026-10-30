"""
Bridge between the supplement dashboard and the unchanged two-stage fractionation
code of this package (fluid_fractionation.py and its engine modules).

The two-stage model reads its inputs from Streamlit session state. This module
prepares those entries from `data/sample_input_data.xlsx` exactly as the project
dashboard did, then calls `render_fluid_fractionation_section()`. Nothing in
the engine modules is modified by this bridge.
"""

import os
import sys

import numpy as np
import pandas as pd
import streamlit as st

P2_DIR = os.path.dirname(os.path.abspath(__file__))                 # this package (engine modules)
ROOT = os.path.dirname(P2_DIR)                                        # repository root
DATA_DIR = os.path.join(ROOT, "data")
INPUT_XLSX = os.path.join(DATA_DIR, "sample_input_data.xlsx")
ARC_XLSX = os.path.join(DATA_DIR, "arc_data.xlsx")

if P2_DIR not in sys.path:
    sys.path.insert(0, P2_DIR)

from CORE_scientific_calculations_engine import load_excel, normalize_headers, process_all_vectorized  # noqa: E402
import arc_data_overlay  # noqa: E402

# arc_data_overlay looks for the arc workbook next to itself; point it at data/ instead
# (os.path.join with an absolute second argument returns that absolute path).
arc_data_overlay.ARC_DATA_FILENAME = ARC_XLSX

SALINITIES = [0.001, 5.0, 10.0, 15.0, 20.0]


def safe_float(value, default):
    try:
        if value is None or (pd.notna(value) is False):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


@st.cache_data
def load_part1_source(_mtime: float):
    """First valid row of the 'Part1' sheet: D_bulk table at the engine's salinity grid and C0."""
    try:
        raw = load_excel(INPUT_XLSX, sheet="Part1")
    except Exception:
        raw = load_excel(INPUT_XLSX, sheet=0)
    df = normalize_headers(raw)
    for pos in range(len(df)):
        row = df.iloc[pos]
        if np.isfinite(safe_float(row.get("T_K"), np.nan)) and np.isfinite(safe_float(row.get("P_bar"), np.nan)):
            results = process_all_vectorized(df.iloc[[pos]], SALINITIES)
            comp = {
                "MO": safe_float(row.get("C0_MO"), 0.025), "Ce": safe_float(row.get("C0_Ce"), 0.772),
                "W": safe_float(row.get("C0_W"), 0.0024), "U": safe_float(row.get("C0_U"), 0.0047),
                "Th": safe_float(row.get("C0_Th"), 0.0137), "Nb": safe_float(row.get("C0_Nb"), 0.21),
                "La": safe_float(row.get("C0_La"), 0.234),
            }
            return results, comp
    return None, None


@st.cache_data
def load_part2_segments(_mtime: float):
    """Segments of the 'Part2_Segments' sheet, each with its own D_bulk table (engine at its T_K)."""
    try:
        raw = load_excel(INPUT_XLSX, sheet="Part2_Segments")
    except Exception:
        return []
    cl = {str(c).strip().lower(): c for c in raw.columns}
    if "rock_id" not in cl and "segment" in cl:
        raw["rock_id"] = raw[cl["segment"]]
    try:
        df = normalize_headers(raw)
    except Exception:
        return []
    # Formula cells (T_K, phi) have no cached value after an openpyxl save; recompute them.
    df["T_K"] = pd.to_numeric(df["T_K"], errors="coerce") if "T_K" in df.columns else np.nan
    if "T_C" in df.columns:
        df["T_K"] = df["T_K"].fillna(pd.to_numeric(df["T_C"], errors="coerce") + 273.15)
    if "phi" not in df.columns:
        df["phi"] = np.nan
    df["phi"] = pd.to_numeric(df["phi"], errors="coerce")
    tk = pd.to_numeric(df["T_K"], errors="coerce")
    df.loc[df["phi"].isna() & tk.notna() & (tk > 0), "phi"] = 10.0 ** (4530.0 * (1.0 / 1273.15 - 1.0 / tk))

    def fin(x):
        try:
            return np.isfinite(float(x))
        except (TypeError, ValueError):
            return False

    segments = []
    for i in range(len(df)):
        r = df.iloc[i]
        name = str(r.get("segment", r.get("rock_id")))
        T_K = safe_float(r.get("T_K"), float("nan"))
        phi_excel = safe_float(r.get("phi"), float("nan"))
        phi_from_excel = np.isfinite(phi_excel)
        if phi_from_excel:
            phi_val = phi_excel
        elif np.isfinite(T_K) and T_K > 0:
            phi_val = 10.0 ** (4530.0 * (1.0 / 1273.15 - 1.0 / T_K))
        else:
            phi_val = float("nan")
        required = [r.get("P_bar"), r.get("dFMQ"), T_K, r.get("mode_cpx"), r.get("mode_grt"), r.get("mode_rut"),
                    r.get("C0_MO"), r.get("C0_Ce")]
        ready = all(fin(x) for x in required) and np.isfinite(phi_val)
        results_df = process_all_vectorized(df.iloc[[i]].copy(), SALINITIES) if ready else None
        C0 = {k: safe_float(r.get(col), None) for k, col in
              {"MO": "C0_MO", "Ce": "C0_Ce", "W": "C0_W", "U": "C0_U", "Th": "C0_Th", "Nb": "C0_Nb", "La": "C0_La"}.items()}
        segments.append({
            "name": name,
            "row": i + 2,
            "T_K": T_K,
            "P_bar": safe_float(r.get("P_bar"), float("nan")),
            "dFMQ": safe_float(r.get("dFMQ"), float("nan")),
            "dFMQ_filter": safe_float(r.get("ΔFMQ_filter"), float("nan")),
            "C0_MO_filter": safe_float(r.get("C0_MO_filter"), float("nan")),
            "mode_rut": safe_float(r.get("mode_rut"), 0.0) * 100.0,
            "phi": phi_val,
            "phi_from_excel": bool(phi_from_excel),
            "ready": bool(ready),
            "results_df": results_df,
            "C0": C0,
            "protolith": str(r.get("rock_id", name)),
            "d98Mo": safe_float(r.get("d98Mo"), None),
            "d98Mo_2SD": safe_float(r.get("d98Mo 2SD"), None),
        })
    return segments


def render_two_stage_model():
    try:
        mtime = os.path.getmtime(INPUT_XLSX)
    except OSError:
        st.error(f"Input workbook not found: {INPUT_XLSX}")
        return
    results_df, comp = load_part1_source(mtime)
    if results_df is None:
        st.error("No valid source row in sheet 'Part1' of data/sample_input_data.xlsx.")
        return
    st.session_state["bali_results_df"] = results_df
    st.session_state["bali_source_params"] = comp
    st.session_state["bali_selected_salinity"] = 10.0
    st.session_state["p2_segments"] = load_part2_segments(mtime)

    st.markdown(
        "Evolution of Mo and δ⁹⁸/⁹⁵Mo in the slab fluid, in two stages. "
        "**Stage 1 — extraction (red):** Rayleigh dehydration of the eclogite segment; the bulk partition "
        "coefficients are those of the Mo-W Thermodynamic model evaluated at the segment's own temperature and ΔFMQ "
        "(with the Zr-in-rutile temperature correction φ), and the fluid–eclogite fractionation factor α is "
        "1.0002–1.0004. The curves follow the cumulative extracted fluid as the fluid fraction F increases. "
        "**Stage 2 — reduction (blue):** the fluid crosses a sediment-buffered redox filter at ΔFMQ_filter, "
        "re-equilibrates with the filter rock and loses Mo; the surviving fluid drifts to heavier δ⁹⁸/⁹⁵Mo "
        "with the reduction fractionation factor α_r. Solid curves use the unmodified N-MORB protolith, dashed "
        "curves the sediment-modified protolith set with the sliders. "
        "Segment conditions (P, T, ΔFMQ, ΔFMQ_filter, modes, C₀, δ⁹⁸/⁹⁵Mo) are read from the sheet "
        "`Part2_Segments` of `data/sample_input_data.xlsx`; lava data from `data/arc_data.xlsx`."
    )

    from fluid_fractionation import render_fluid_fractionation_section
    render_fluid_fractionation_section()
