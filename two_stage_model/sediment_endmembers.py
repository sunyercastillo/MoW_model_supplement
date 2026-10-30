"""
Sediment end-member compositions, mixing lines and the end-member markers drawn
on the δ⁹⁸/⁹⁵Mo plots of the two-stage fractionation model.

Contents: carbonate and hemipelagic sediment end-members with δ⁹⁸/⁹⁵Mo,
carbonate–hemipelagic binary mixing (element ratios and isotope hyperbola),
sediment-mix–depleted-mantle mixing, and Plotly overlay helpers.
"""

import numpy as np
import plotly.graph_objects as go

from .partitioning_engine import BaliConstants as _BaliConstants

_bc = _BaliConstants()
# Depleted mantle and N-MORB (µg/g), Table B1 of the article (shared with the partitioning engine)
_KEYS = {"Mo": "MO", "Ce": "Ce", "W": "W", "U": "U", "Th": "Th", "Nb": "Nb", "La": "La"}
DM = {k: _bc.dm_concentrations[v] for k, v in _KEYS.items()}
NMORB = {k: _bc.morb_concentrations[v] for k, v in _KEYS.items()}

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
# δ⁹⁸/⁹⁵Mo values (‰, NIST 3134)
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

# tick marks on mixing lines, smaller than the data markers
_TICK_SIZE = 6

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
#   δ⁹⁸/⁹⁵Mo_mix = (R_A · [Mo]_A · f_A + R_B · [Mo]_B · f_B) / [Mo]_mix
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
# Ratio keys (lower-case aliases accepted)
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
# Overlays for δ⁹⁸/⁹⁵Mo-vs-ratio plots
# ══════════════════════════════════════════════════════════════════════════════

def add_sediment_overlay_xonly(fig, x_ratio_key):
    """Add all sediment overlays to a δ⁹⁸/⁹⁵Mo-vs-ratio plot."""
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
        marker=dict(size=_TICK_SIZE, color=COLOR_CARBONATES, symbol="circle",
                    line=dict(color="black", width=0.8)),
        text=_ch_labels2,
        textposition="top center",
        textfont=dict(size=7, color=COLOR_CARBONATES),
        showlegend=False,
        hovertemplate="f_hemi = %{text}<extra></extra>",
    ))

    # ── Endmember stars at (ratio, δ⁹⁸/⁹⁵Mo) ──
    fig.add_trace(go.Scatter(
        x=[CARB_RATIOS[xk]], y=[D98MO_CARBONATES],
        mode="markers",
        marker=dict(size=MARKER_SIZE, color=COLOR_CARBONATES,
                    symbol=MARKER_SYMBOL, line=dict(color="black", width=0.8)),
        name="Subd. Carbonates", legendgroup="sed_carb", showlegend=False,
        hovertemplate=(
            f"<b>Subd. Carbonates</b><br>"
            f"{xk}: %{{x:.3f}}<br>δ⁹⁸/⁹⁵Mo: %{{y:.2f}}‰<extra></extra>"
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
            f"{xk}: %{{x:.3f}}<br>δ⁹⁸/⁹⁵Mo: %{{y:.2f}}‰<extra></extra>"
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
        marker=dict(size=_TICK_SIZE, color=COLOR_TICK, symbol="circle",
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
            f"{xk}: %{{x:.3f}}<br>δ⁹⁸/⁹⁵Mo: %{{y:.2f}}‰<extra></extra>"
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
            f"{xk}: %{{x:.3f}}<br>δ⁹⁸/⁹⁵Mo: %{{y:.2f}}‰<extra></extra>"
        ),
    ))

    return fig


def add_endmember_markers_xonly(fig, x_ratio_key):
    """Add ONLY endmember markers (no mixing lines) to a δ⁹⁸/⁹⁵Mo-vs-ratio plot."""
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
            f"{xk}: %{{x:.3f}}<br>δ⁹⁸/⁹⁵Mo: %{{y:.2f}}‰<extra></extra>"
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
            f"{xk}: %{{x:.3f}}<br>δ⁹⁸/⁹⁵Mo: %{{y:.2f}}‰<extra></extra>"
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
            f"{xk}: %{{x:.3f}}<br>δ⁹⁸/⁹⁵Mo: %{{y:.2f}}‰<extra></extra>"
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
            f"{xk}: %{{x:.3f}}<br>δ⁹⁸/⁹⁵Mo: %{{y:.2f}}‰<extra></extra>"
        ),
    ))

    return fig
