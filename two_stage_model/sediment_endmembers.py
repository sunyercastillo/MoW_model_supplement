"""
Sediment endmember compositions, mixing-line calculations, and Plotly overlays.

Contents
--------
• Two subducted-sediment endmembers (Carbonates, Hemipelagic) with δ⁹⁸Mo
• Carb ↔ Hemi binary mixing (element ratios + δ⁹⁸Mo isotope hyperbola)
• 50 % Sed-mix ↔ DM mixing (element ratios + δ⁹⁸Mo isotope hyperbola)
• Public overlay functions for Part 1 (ratio-ratio) and Part 2 (δ⁹⁸Mo-ratio) plots

This module does NOT modify any core scientific calculation files.
"""

import numpy as np
import plotly.graph_objects as go

# ══════════════════════════════════════════════════════════════════════════════
# Depleted Mantle (DM) composition — read from the core engine (not hardcoded)
# ══════════════════════════════════════════════════════════════════════════════

from CORE_scientific_calculations_engine import BaliConstants as _BaliConstants

_bc = _BaliConstants()
_DM_RAW = _bc.dm_concentrations              # keys: "MO", "Ce", "W", "U", "Th", "Nb", "La"

DM = {
    "Mo": _DM_RAW["MO"],    # 0.025  µg/g — Salter & Stracke (2004)
    "Ce": _DM_RAW["Ce"],    # 0.772  µg/g — Salter & Stracke (2004)
    "W":  _DM_RAW["W"],     # 0.0024 µg/g — König et al. (2011)
    "U":  _DM_RAW["U"],     # 0.0047 µg/g — Salter & Stracke (2004)
    "Th": _DM_RAW["Th"],    # 0.0137 µg/g — Salter & Stracke (2004)
    "Nb": _DM_RAW["Nb"],    # 0.21   µg/g — Salter & Stracke (2004)
    "La": _DM_RAW["La"],    # 0.234  µg/g — Salter & Stracke (2004)
}

_MORB_RAW = _bc.morb_concentrations

NMORB = {
    "Mo": _MORB_RAW["MO"],   # 0.19  µg/g — Bezard et al. (2016)
    "Ce": _MORB_RAW["Ce"],   # 5.6   µg/g — Bezard et al. (2016)
    "W":  _MORB_RAW["W"],    # 0.040 µg/g — König et al. (2011)
    "U":  _MORB_RAW["U"],    # 0.047 µg/g — Gale et al. (2013)
    "Th": _MORB_RAW["Th"],   # 0.087 µg/g — Gale et al. (2013)
    "Nb": _MORB_RAW["Nb"],   # 2.33  µg/g — Gale et al. (2013)
    "La": _MORB_RAW["La"],   # 2.50  µg/g — Gale et al. (2013)
}

# ══════════════════════════════════════════════════════════════════════════════
# Endmember compositions (ppm = µg/g)
# ══════════════════════════════════════════════════════════════════════════════

# Arithmetic mean of three samples per lithological unit, ODP Leg 170 Site 1039 (EA-1); carbonate excludes the calcareous breccia at 367 mbsf.
SUBDUCTED_CARBONATES = {
    "Nb": 1.317,  "Ba": 1796.654, "La": 5.254,  "Ce": 2.960,
    "Th": 0.222,  "Nd": 3.848,    "Yb": 0.748,  "Sr": 1338.700,
    "Hf": 0.230,  "Pb": 1.670,    "U":  0.684,  "Mo": 1.087,
    "W":  0.187,  "Ta": 0.077,
}

# Arithmetic mean of three samples per lithological unit, ODP Leg 170 Site 1039 (EA-1); carbonate excludes the calcareous breccia at 367 mbsf.
SUBDUCTED_HEMIPELAGIC = {
    "Nb": 5.910,  "Ba": 1710.165, "La": 12.862, "Ce": 25.774,
    "Th": 2.656,  "Nd": 13.635,   "Yb": 2.031,  "Sr": 236.181,
    "Hf": 2.488,  "Pb": 25.942,   "U":  4.117,  "Mo": 5.626,
    "W":  0.874,  "Ta": 0.360,
}

# ══════════════════════════════════════════════════════════════════════════════
# δ⁹⁸Mo values (‰, NIST 3134)
# ══════════════════════════════════════════════════════════════════════════════

D98MO_CARBONATES  = +0.87
D98MO_HEMIPELAGIC = -2.52
D98MO_DM          = -0.21
D98MO_NMORB       = -0.21

# ══════════════════════════════════════════════════════════════════════════════
# Colours, markers, tick-mark sizes
# ══════════════════════════════════════════════════════════════════════════════

COLOR_CARBONATES  = "#87CEEB"   # light blue
COLOR_HEMIPELAGIC = "#8B4513"   # brown
COLOR_MORB_SED    = "#AAAAAA"   # grey for MORB–Sed mixing
COLOR_DM_STAR     = "#FFD700"   # yellow for DM endmember star
COLOR_NMORB       = "#AAAAAA"   # grey for N-MORB circle
COLOR_TICK        = "#888888"   # grey tick marks
MARKER_SYMBOL     = "star"
MARKER_SIZE       = 12

# 33 % smaller than normal data-point markers in each plot section
_TICK_SIZE_PART1 = 4     # Part 1 normal marker = 6  → 6 × 0.67 ≈ 4
_TICK_SIZE_PART2 = 6     # Part 2 normal marker = 9  → 9 × 0.67 ≈ 6

# ══════════════════════════════════════════════════════════════════════════════
# Ratio helpers
# ══════════════════════════════════════════════════════════════════════════════

_RATIO_ELEMENTS = {
    "Mo_Ce": ("Mo", "Ce"),
    "Ce_Mo": ("Ce", "Mo"),
    "U_Th":  ("U",  "Th"),
    "W_Th":  ("W",  "Th"),
    "Mo_W":  ("Mo", "W"),
    "Nb_La": ("Nb", "La"),
    "W_Ce":  ("W",  "Ce"),
}

_ELEMENTS_FOR_RATIOS = ("Mo", "Ce", "U", "Th", "W", "Nb", "La")


def _ratio(comp, num_key, den_key):
    n = comp.get(num_key, 0.0)
    d = comp.get(den_key, 0.0)
    return n / d if d != 0 else float("nan")


def _endmember_ratios(comp):
    return {k: _ratio(comp, n, d) for k, (n, d) in _RATIO_ELEMENTS.items()}


def _safe_array_ratio(a, b):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(b != 0, a / b, np.nan)


CARB_RATIOS  = _endmember_ratios(SUBDUCTED_CARBONATES)
HEMI_RATIOS  = _endmember_ratios(SUBDUCTED_HEMIPELAGIC)
DM_RATIOS    = _endmember_ratios(DM)
NMORB_RATIOS = _endmember_ratios(NMORB)

# ══════════════════════════════════════════════════════════════════════════════
# Carb–Hemi binary mixing  (f_hemi = 0 → pure Carbonates, 1 → pure Hemipelagic)
# ══════════════════════════════════════════════════════════════════════════════

_F_HEMI = np.linspace(0, 1, 101)


def _mix_two_arrays(comp_a, comp_b, f_b, elements):
    """Linear element mixing: mix = comp_a·(1-f_b) + comp_b·f_b."""
    return {
        el: comp_a.get(el, 0.0) * (1 - f_b) + comp_b.get(el, 0.0) * f_b
        for el in elements
    }


_MIX_CH = _mix_two_arrays(SUBDUCTED_CARBONATES, SUBDUCTED_HEMIPELAGIC,
                           _F_HEMI, _ELEMENTS_FOR_RATIOS)

MIXING_RATIOS_CH = {
    k: _safe_array_ratio(_MIX_CH[n], _MIX_CH[d])
    for k, (n, d) in _RATIO_ELEMENTS.items()
}

# Isotope mixing (hyperbolic):
#   δ⁹⁸Mo_mix = (R_A · [Mo]_A · f_A + R_B · [Mo]_B · f_B) / [Mo]_mix
MIXING_D98MO_CH = (
    D98MO_CARBONATES  * SUBDUCTED_CARBONATES["Mo"]  * (1 - _F_HEMI) +
    D98MO_HEMIPELAGIC * SUBDUCTED_HEMIPELAGIC["Mo"] * _F_HEMI
) / _MIX_CH["Mo"]

_MASK_CARB = _F_HEMI <= 0.5
_MASK_HEMI = _F_HEMI >= 0.5      # overlap at f=0.5 ensures continuity

# Carb-Hemi tick marks: 25 %, 50 %, 75 % only
_TICK_IDX_CH = [25, 50, 75]

# Backward-compat alias used by external callers
MIXING_RATIOS = MIXING_RATIOS_CH

# ══════════════════════════════════════════════════════════════════════════════
# 50 % Sediment Mix (equal parts Carbonate + Hemipelagic)
# ══════════════════════════════════════════════════════════════════════════════

SED_50MIX = {
    el: 0.5 * SUBDUCTED_CARBONATES.get(el, 0.0)
      + 0.5 * SUBDUCTED_HEMIPELAGIC.get(el, 0.0)
    for el in _ELEMENTS_FOR_RATIOS
}
SED_50MIX_RATIOS = {
    k: _ratio(SED_50MIX, n, d) for k, (n, d) in _RATIO_ELEMENTS.items()
}

D98MO_50MIX = (
    D98MO_CARBONATES  * SUBDUCTED_CARBONATES["Mo"]  * 0.5 +
    D98MO_HEMIPELAGIC * SUBDUCTED_HEMIPELAGIC["Mo"] * 0.5
) / SED_50MIX["Mo"]


def make_mixed_sediment(f_hemi):
    """Return Carb–Hemi blend at f_hemi (0 = pure Carbonates, 1 = pure Hemipelagic)."""
    return _mix_two_arrays(SUBDUCTED_CARBONATES, SUBDUCTED_HEMIPELAGIC,
                           f_hemi, _ELEMENTS_FOR_RATIOS)


# ══════════════════════════════════════════════════════════════════════════════
# DM–Sediment mixing  (f_dm = 0 → pure 50 % sed, 1 → pure DM)
# ══════════════════════════════════════════════════════════════════════════════

_F_DM = np.linspace(0, 1, 101)

_MIX_MS = _mix_two_arrays(SED_50MIX, DM, _F_DM, _ELEMENTS_FOR_RATIOS)

MIXING_RATIOS_MS = {
    k: _safe_array_ratio(_MIX_MS[n], _MIX_MS[d])
    for k, (n, d) in _RATIO_ELEMENTS.items()
}

# Isotope mixing (hyperbolic)
MIXING_D98MO_MS = (
    D98MO_50MIX * SED_50MIX["Mo"] * (1 - _F_DM) +
    D98MO_DM    * DM["Mo"]        * _F_DM
) / _MIX_MS["Mo"]

# DM-Sed tick marks: every 10 %
_TICK_IDX_MS = [90, 98]                 # f_DM = 90 % and 98 %

# ══════════════════════════════════════════════════════════════════════════════
# Key mapping  (dashboard Part 1 uses "u_th", Part 2 uses "U_Th")
# ══════════════════════════════════════════════════════════════════════════════

_RATIO_KEY_MAP = {
    "u_th": "U_Th",  "mo_ce": "Mo_Ce",  "ce_mo": "Ce_Mo",  "w_th": "W_Th",
    "mo_w": "Mo_W",  "nb_la": "Nb_La",  "w_ce": "W_Ce",
    "U_Th": "U_Th",  "Mo_Ce": "Mo_Ce",  "Ce_Mo": "Ce_Mo",  "W_Th": "W_Th",
    "Mo_W": "Mo_W",  "Nb_La": "Nb_La",  "W_Ce": "W_Ce",
}


def _resolve(key):
    return _RATIO_KEY_MAP.get(key, key)


# ══════════════════════════════════════════════════════════════════════════════
# PUBLIC API — Part 1 (both axes are element ratios)
# ══════════════════════════════════════════════════════════════════════════════

def add_sediment_overlay(fig, x_ratio_key, y_ratio_key):
    """Add endmember overlays (DM star, N-MORB circle) to a Part 1 ratio-vs-ratio plot.

    The Subd. Carbonates / Subd. Hemipelagic stars are deliberately NOT drawn in
    Part 1 (hidden per design decision); Part 2 keeps them via the _xonly variants.
    """
    xk = _resolve(x_ratio_key)
    yk = _resolve(y_ratio_key)

    if xk not in MIXING_RATIOS_CH or yk not in MIXING_RATIOS_CH:
        return fig

    # ── DM endmember star (yellow) ──
    fig.add_trace(go.Scatter(
        x=[DM_RATIOS[xk]], y=[DM_RATIOS[yk]],
        mode="markers",
        marker=dict(size=MARKER_SIZE, color=COLOR_DM_STAR,
                    symbol=MARKER_SYMBOL, line=dict(color="black", width=0.8)),
        name="DM", legendgroup="dm", showlegend=False,
        hovertemplate=(
            f"<b>Depleted Mantle</b><br>"
            f"{xk}: %{{x:.3f}}<br>{yk}: %{{y:.3f}}<extra></extra>"
        ),
    ))

    # ── N-MORB crosshair dashed lines (behind all other traces) ──
    fig.add_vline(x=NMORB_RATIOS[xk], line_dash='dash', line_width=1.15,
                  line_color='grey', layer='below')
    fig.add_hline(y=NMORB_RATIOS[yk], line_dash='dash', line_width=1.15,
                  line_color='grey', layer='below')

    return fig


# ══════════════════════════════════════════════════════════════════════════════
# PUBLIC API — Part 2 (x = element ratio, y = δ⁹⁸Mo)
# ══════════════════════════════════════════════════════════════════════════════

def add_sediment_overlay_xonly(fig, x_ratio_key):
    """Add all sediment overlays to a Part 2 δ⁹⁸Mo-vs-ratio plot."""
    xk = _resolve(x_ratio_key)
    if xk not in MIXING_RATIOS_CH:
        return fig

    # ── Carb–Hemi isotope mixing line (hyperbolic, two-colour, solid) ──
    x_ch = MIXING_RATIOS_CH[xk]
    y_ch = MIXING_D98MO_CH

    fig.add_trace(go.Scatter(
        x=x_ch[_MASK_CARB], y=y_ch[_MASK_CARB],
        mode="lines",
        line=dict(color=COLOR_CARBONATES, width=1.5),
        showlegend=False, hoverinfo="skip",
    ))
    fig.add_trace(go.Scatter(
        x=x_ch[_MASK_HEMI], y=y_ch[_MASK_HEMI],
        mode="lines",
        line=dict(color=COLOR_HEMIPELAGIC, width=1.5),
        showlegend=False, hoverinfo="skip",
    ))

    # ── Carb–Hemi tick marks (25 / 50 / 75 %, light blue, 33 % smaller) ──
    _ch_labels2 = [f"{int(_F_HEMI[i]*100)}%" for i in _TICK_IDX_CH]
    fig.add_trace(go.Scatter(
        x=[x_ch[i] for i in _TICK_IDX_CH],
        y=[y_ch[i] for i in _TICK_IDX_CH],
        mode="markers+text",
        marker=dict(size=_TICK_SIZE_PART2, color=COLOR_CARBONATES, symbol="circle",
                    line=dict(color="black", width=0.8)),
        text=_ch_labels2,
        textposition="top center",
        textfont=dict(size=7, color=COLOR_CARBONATES),
        showlegend=False,
        hovertemplate="f_hemi = %{text}<extra></extra>",
    ))

    # ── Endmember stars at (ratio, δ⁹⁸Mo) ──
    fig.add_trace(go.Scatter(
        x=[CARB_RATIOS[xk]], y=[D98MO_CARBONATES],
        mode="markers",
        marker=dict(size=MARKER_SIZE, color=COLOR_CARBONATES,
                    symbol=MARKER_SYMBOL, line=dict(color="black", width=0.8)),
        name="Subd. Carbonates", legendgroup="sed_carb", showlegend=False,
        hovertemplate=(
            f"<b>Subd. Carbonates</b><br>"
            f"{xk}: %{{x:.3f}}<br>δ⁹⁸Mo: %{{y:.2f}}‰<extra></extra>"
        ),
    ))
    fig.add_trace(go.Scatter(
        x=[HEMI_RATIOS[xk]], y=[D98MO_HEMIPELAGIC],
        mode="markers",
        marker=dict(size=MARKER_SIZE, color=COLOR_HEMIPELAGIC,
                    symbol=MARKER_SYMBOL, line=dict(color="black", width=0.8)),
        name="Subd. Hemipelagic", legendgroup="sed_hemi", showlegend=False,
        hovertemplate=(
            f"<b>Subd. Hemipelagic</b><br>"
            f"{xk}: %{{x:.3f}}<br>δ⁹⁸Mo: %{{y:.2f}}‰<extra></extra>"
        ),
    ))

    # ── MORB–Sed isotope mixing line (solid grey, hyperbolic) ──
    x_ms = MIXING_RATIOS_MS[xk]
    y_ms = MIXING_D98MO_MS

    fig.add_trace(go.Scatter(
        x=x_ms, y=y_ms,
        mode="lines",
        line=dict(color=COLOR_MORB_SED, width=1.5),
        showlegend=False, hoverinfo="skip",
    ))

    # ── DM–Sed tick marks (every 10 %, grey, 33 % smaller) ──
    _ms_labels2 = [f"{int(_F_DM[i]*100)}%" for i in _TICK_IDX_MS]
    fig.add_trace(go.Scatter(
        x=x_ms[_TICK_IDX_MS], y=y_ms[_TICK_IDX_MS],
        mode="markers+text",
        marker=dict(size=_TICK_SIZE_PART2, color=COLOR_TICK, symbol="circle",
                    line=dict(color="black", width=0.8)),
        text=_ms_labels2,
        textposition="top center",
        textfont=dict(size=7, color=COLOR_TICK),
        showlegend=False,
        hovertemplate="f_DM = %{text}<extra></extra>",
    ))

    # ── DM endmember star (yellow) ──
    fig.add_trace(go.Scatter(
        x=[DM_RATIOS[xk]], y=[D98MO_DM],
        mode="markers",
        marker=dict(size=MARKER_SIZE, color=COLOR_DM_STAR,
                    symbol=MARKER_SYMBOL, line=dict(color="black", width=0.8)),
        name="DM", legendgroup="dm", showlegend=False,
        hovertemplate=(
            f"<b>Depleted Mantle</b><br>"
            f"{xk}: %{{x:.3f}}<br>δ⁹⁸Mo: %{{y:.2f}}‰<extra></extra>"
        ),
    ))

    # ── N-MORB circle (grey, black outline) ──
    fig.add_trace(go.Scatter(
        x=[NMORB_RATIOS[xk]], y=[D98MO_NMORB],
        mode="markers",
        marker=dict(size=MARKER_SIZE, color=COLOR_NMORB,
                    symbol="circle", line=dict(color="black", width=0.8)),
        name="N-MORB", legendgroup="nmorb", showlegend=False,
        hovertemplate=(
            f"<b>N-MORB</b><br>"
            f"{xk}: %{{x:.3f}}<br>δ⁹⁸Mo: %{{y:.2f}}‰<extra></extra>"
        ),
    ))

    return fig


def add_endmember_markers_xonly(fig, x_ratio_key):
    """Add ONLY endmember markers (no mixing lines) to a Part 2 δ⁹⁸Mo-vs-ratio plot."""
    xk = _resolve(x_ratio_key)
    if xk not in MIXING_RATIOS_CH:
        return fig

    # Carb star
    fig.add_trace(go.Scatter(
        x=[CARB_RATIOS[xk]], y=[D98MO_CARBONATES],
        mode="markers",
        marker=dict(size=MARKER_SIZE, color=COLOR_CARBONATES,
                    symbol=MARKER_SYMBOL, line=dict(color="black", width=0.8)),
        name="Subd. Carbonates", legendgroup="sed_carb", showlegend=False,
        hovertemplate=(
            f"<b>Subd. Carbonates</b><br>"
            f"{xk}: %{{x:.3f}}<br>δ⁹⁸Mo: %{{y:.2f}}‰<extra></extra>"
        ),
    ))
    # Hemi star
    fig.add_trace(go.Scatter(
        x=[HEMI_RATIOS[xk]], y=[D98MO_HEMIPELAGIC],
        mode="markers",
        marker=dict(size=MARKER_SIZE, color=COLOR_HEMIPELAGIC,
                    symbol=MARKER_SYMBOL, line=dict(color="black", width=0.8)),
        name="Subd. Hemipelagic", legendgroup="sed_hemi", showlegend=False,
        hovertemplate=(
            f"<b>Subd. Hemipelagic</b><br>"
            f"{xk}: %{{x:.3f}}<br>δ⁹⁸Mo: %{{y:.2f}}‰<extra></extra>"
        ),
    ))
    # DM star
    fig.add_trace(go.Scatter(
        x=[DM_RATIOS[xk]], y=[D98MO_DM],
        mode="markers",
        marker=dict(size=MARKER_SIZE, color=COLOR_DM_STAR,
                    symbol=MARKER_SYMBOL, line=dict(color="black", width=0.8)),
        name="DM", legendgroup="dm", showlegend=False,
        hovertemplate=(
            f"<b>Depleted Mantle</b><br>"
            f"{xk}: %{{x:.3f}}<br>δ⁹⁸Mo: %{{y:.2f}}‰<extra></extra>"
        ),
    ))
    # N-MORB circle
    fig.add_trace(go.Scatter(
        x=[NMORB_RATIOS[xk]], y=[D98MO_NMORB],
        mode="markers",
        marker=dict(size=MARKER_SIZE, color=COLOR_NMORB,
                    symbol="circle", line=dict(color="black", width=0.8)),
        name="N-MORB", legendgroup="nmorb", showlegend=False,
        hovertemplate=(
            f"<b>N-MORB</b><br>"
            f"{xk}: %{{x:.3f}}<br>δ⁹⁸Mo: %{{y:.2f}}‰<extra></extra>"
        ),
    ))

    return fig


# ══════════════════════════════════════════════════════════════════════════════
# PUBLIC API — Slab sediment influence (Part 1 ratio-ratio plots)
# ══════════════════════════════════════════════════════════════════════════════

# Fluid fraction grid (same as existing wedge mixing in the dashboard)
_F_FLUID_SLAB = [0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.9]


def compute_slab_fluid_endmember(f_sed, sed_comp, morb_comp, D_bulk_dict):
    """Return fluid concentrations for a slab = (1-f_sed)·MORB + f_sed·sediment.

    Parameters
    ----------
    f_sed : float
        Fraction of sediment in the slab source (0 = pure MORB, 1 = pure sed).
    sed_comp : dict
        Sediment endmember concentrations (ppm).
    morb_comp : dict
        N-MORB concentrations (ppm).
    D_bulk_dict : dict
        Bulk partition coefficients keyed by element name
        (e.g. {"Mo": 0.05, "Ce": 2.1, ...}).

    Returns
    -------
    dict  Element → fluid concentration (ppm).
    """
    elements = ("Mo", "Ce", "W", "U", "Th", "Nb", "La")
    C_fluid = {}
    for el in elements:
        C0_mixed = (1 - f_sed) * morb_comp.get(el, 0.0) + f_sed * sed_comp.get(el, 0.0)
        D = D_bulk_dict.get(el, 0.0)
        C_fluid[el] = C0_mixed / D if D > 0 else float("nan")
    return C_fluid


def add_slab_influence_mixing_line(
    fig, x_ratio_key, y_ratio_key,
    D_bulk_dict, sed_comp, f_sed, dm_comp,
    salinity_color, morb_comp=None,
):
    """Add a slab-influence mixing line to a Part 1 ratio-ratio plot.

    The line is the wedge-mixing curve (fluid ↔ DM) computed from the
    slab fluid endmember, which is a blend of a MORB protolith and a sediment
    component. Style: dashed, same colour as the reference salinity.

    Parameters
    ----------
    fig : go.Figure
    x_ratio_key : str   e.g. "u_th"
    y_ratio_key : str   e.g. "mo_ce"
    D_bulk_dict : dict  Bulk D values from excel_results_df row (element keys)
    sed_comp : dict     Sediment endmember (one of SUBDUCTED_CARBONATES etc.)
    f_sed : float       Sediment fraction in slab (0–1)
    dm_comp : dict      Depleted Mantle composition (DM)
    salinity_color : str  Hex colour string for line and markers
    morb_comp : dict or None
        MORB protolith concentrations. Defaults to canonical N-MORB if None.

    Returns
    -------
    go.Figure
    """
    xk = _resolve(x_ratio_key)
    yk = _resolve(y_ratio_key)
    if xk not in _RATIO_ELEMENTS or yk not in _RATIO_ELEMENTS:
        return fig

    x_num, x_den = _RATIO_ELEMENTS[xk]
    y_num, y_den = _RATIO_ELEMENTS[yk]

    # Fluid endmember from slab (mixed source)
    _morb = morb_comp if morb_comp is not None else NMORB
    C_fluid = compute_slab_fluid_endmember(f_sed, sed_comp, _morb, D_bulk_dict)

    # Binary mixing: C_mix = f·C_fluid + (1-f)·C_DM
    x_vals, y_vals, f_labels = [], [], []
    for f in _F_FLUID_SLAB:
        C_mix = {
            el: f * C_fluid.get(el, 0.0) + (1 - f) * dm_comp.get(el, 0.0)
            for el in (x_num, x_den, y_num, y_den)
        }
        xv = C_mix[x_num] / C_mix[x_den] if C_mix[x_den] != 0 else float("nan")
        yv = C_mix[y_num] / C_mix[y_den] if C_mix[y_den] != 0 else float("nan")
        x_vals.append(xv)
        y_vals.append(yv)
        f_labels.append(f"{f*100:.1f}%")

    sed_pct = int(round(f_sed * 100))

    # Main dashed line (slab influence)
    fig.add_trace(go.Scatter(
        x=x_vals, y=y_vals,
        mode="lines",
        line=dict(color=salinity_color, width=1.5, dash='dash'),
        showlegend=False,
        hovertemplate=(
            f"<b>Slab influence ({sed_pct}% sed)</b><br>"
            f"f_fluid = %{{text}}<br>"
            f"{xk}: %{{x:.3f}}<br>{yk}: %{{y:.3f}}<extra></extra>"
        ),
        text=f_labels,
    ))

    # Tick marker + "1%" label at f_fluid = 0.01
    idx_1pct = _F_FLUID_SLAB.index(0.01)
    fig.add_trace(go.Scatter(
        x=[x_vals[idx_1pct]], y=[y_vals[idx_1pct]],
        mode="markers+text",
        marker=dict(size=6, symbol="circle", color=salinity_color,
                    line=dict(color='black', width=0.5)),
        text=["1%"],
        textposition="top center",
        textfont=dict(color=salinity_color, size=10,
                      family="Arial, sans-serif", weight="bold"),
        showlegend=False,
        hoverinfo="skip",
    ))

    # Label at end of line (last point)
    fig.add_trace(go.Scatter(
        x=[x_vals[-1]], y=[y_vals[-1]],
        mode="text",
        text=[f"{sed_pct}%sed"],
        textposition="top center",
        textfont=dict(color=salinity_color, size=8, family="Arial, sans-serif"),
        showlegend=False,
        hoverinfo="skip",
    ))

    return fig
