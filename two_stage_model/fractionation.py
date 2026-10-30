"""
Two-stage fractionation model.

Stage 1 — Rayleigh dehydration of an eclogite slab segment: Mo, W, U, Th, Ce,
Nb, La and δ⁹⁸/⁹⁵Mo of the cumulative extracted fluid as a function of the fluid
fraction F, with bulk partition coefficients from the partitioning engine at the
segment's own temperature and ΔFMQ and a fluid–eclogite fractionation factor α.

Stage 2 — redox filter: the fluid re-equilibrates with a sediment buffered at
ΔFMQ_filter, loses Mo and drifts to heavier δ⁹⁸/⁹⁵Mo (fractionation factor α_r).

Inputs are prepared by bridge.py from the workbooks in data/; the public entry
point is render_fluid_fractionation_section().
"""

import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from typing import Dict, List, Tuple, Optional

from .lava_data import add_arc_overlay_to_figure, load_arc_data_from_file
from .sediment_endmembers import (
    make_mixed_sediment,
    D98MO_CARBONATES, D98MO_HEMIPELAGIC,
    SUBDUCTED_CARBONATES, SUBDUCTED_HEMIPELAGIC,
)

# ============================================================================
# CONSTANTS
# ============================================================================
# φ (Zr-in-rutile temperature correction) is read per segment from the
# 'Part2_Segments' sheet (column `phi`).
W_D_RUT_FIXED = 1.250       # Fixed W rutile partition coefficient (NOT corrected by φ)

# ============================================================================
# F grid
# ============================================================================
F_LIST_DEFAULT = [
    0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
    0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19,
    0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29,
    0.30, 0.33, 0.35, 0.38, 0.40, 0.43, 0.45, 0.48, 0.50,
    0.60, 0.70, 0.80, 0.90, 0.99
]

# Highlight F values for labeling (only 0.1, 0.2, 0.3)
F_HIGHLIGHT = [0.10, 0.20, 0.30]


# ============================================================================
# STAGE-2 REDOX-FILTER TRANSPORT CONSTANTS (fluid LOSES Mo; removal EMERGENT)
# ============================================================================
# Stage 2 is a redox filter: the slab-derived fluid percolates through a
# reactive sediment horizon buffered at an ABSOLUTE fO₂ (ΔFMQ_filter — a
# petrological buffer floor, e.g. graphite/CCO-class metasediment ≈ FMQ
# −0.5…−4). Mechanism = fluid–rock RE-EQUILIBRATION (partitioning): the fluid
# exchanges Mo with the filter rock and its exit value is rock-buffered,
#     C_eq = C_rock^filter / D_bulk(T, ΔFMQ_filter)      (clipped ≤ C_start)
# valid for the dilute Stage-1 fluid (×600–1000 below the Bali solubility
# ceiling, so a saturation-limited law does not apply). The exchange
# progress f_ex (0→1) spans no-exchange → full equilibration; the clip
# guarantees the sediment never ADDS net Mo. Note C_exit is NOT zero: a
# Mo-bearing rock buffers the fluid at C_rock/D — equilibrium forbids total
# stripping (and the nonzero floor is what the lavas record).
# It branches from the Stage-1 CUMULATIVE extracted fluid at F = F_BRANCH.
F_BRANCH = 0.20                     # Stage-1 F at which Stage-2 reduction branches
# (coincides with the F=0.20 highlighted marker, F_HIGHLIGHT)
ALPHA_EXT_REF = 1.0003              # extraction α used to fix the Stage-2 start point
# α_r family: Stage-1 values + a near-zero option (1.00002, "barely fractionating"
# reduction — the value class required for the same-slab/Nicaragua-filter scenario
# to survive Rayleigh removal). A custom α_r can be added from the dashboard.
ALPHA_R_DEFAULT = (1.00002, 1.0002, 1.0003, 1.0004)
ALPHA_R_REF = 1.0003               # default-displayed reduction α_r
EXP_MO, EXP_W, EXP_U = 0.44, 0.07, 0.1433    # Bali fO₂ exponents (do not change)
# Removal markers: 25/50/75/90 % Mo removed ⇒ g = 0.75/0.50/0.25/0.10
# (rows at these exact survivals are inserted analytically when reachable)
G_HIGHLIGHT = [0.75, 0.50, 0.25, 0.10]


# ============================================================================
# CORE COMPUTATION FUNCTIONS
# ============================================================================

@st.cache_data(show_spinner=False)
def compute_fluid_fractionation(
    C0_MO: float,
    C0_Ce: float,
    MO_Dbulk: float,
    CE_Dbulk: float,
    C0_U: float = None,
    C0_Th: float = None,
    U_Dbulk: float = None,
    TH_Dbulk: float = None,
    C0_W: float = None,
    W_Dbulk: float = None,
    C0_Nb: float = None,
    C0_La: float = None,
    NB_Dbulk: float = None,
    LA_Dbulk: float = None,
    delta0: float = -0.21,
    alphas: Tuple[float, ...] = (1.0002, 1.0003, 1.0004),
    F_list: List[float] = None
) -> Dict[float, pd.DataFrame]:
    """
    Compute cumulative Rayleigh fluid fractionation model for Mo isotopes, Mo/Ce, U/Th, W/Th, and Mo/W.
    
    Parameters
    ----------
    C0_MO : float
        Initial Mo concentration in source rock (ppm)
    C0_Ce : float
        Initial Ce concentration in source rock (ppm)
    MO_Dbulk : float
        Bulk solid/fluid partition coefficient for Mo
    CE_Dbulk : float
        Bulk solid/fluid partition coefficient for Ce
    C0_U : float, optional
        Initial U concentration in source rock (ppm)
    C0_Th : float, optional
        Initial Th concentration in source rock (ppm)
    U_Dbulk : float, optional
        Bulk solid/fluid partition coefficient for U
    TH_Dbulk : float, optional
        Bulk solid/fluid partition coefficient for Th
    C0_W : float, optional
        Initial W concentration in source rock (ppm)
    W_Dbulk : float, optional
        Bulk solid/fluid partition coefficient for W
    delta0 : float, optional
        Initial δ⁹⁸/⁹⁵Mo isotope value (‰), default -0.21
    alphas : tuple of float, optional
        Fractionation factors to compute, default (1.0002, 1.0003, 1.0004)
    F_list : list of float, optional
        Cumulative mass fractions of fluid extracted, default F_LIST_DEFAULT
    
    Returns
    -------
    dict
        Dictionary mapping alpha -> pandas.DataFrame with columns:
        F, cs_Mo, cs_Ce, cs_U, cs_Th, cs_W, fs, Rs, delta_s, 
        cf_cum_Mo, cf_cum_Ce, cf_cum_U, cf_cum_Th, cf_cum_W,
        MoCe_f_cum, UTh_f_cum, WTh_f_cum, MoW_f_cum, Rf_cum, delta_f_cum
    
    Raises
    ------
    ValueError
        If input parameters are invalid or validation fails
    """
    if F_list is None:
        F_list = F_LIST_DEFAULT.copy()
    
    # Validate inputs for Mo and Ce (required)
    if C0_MO <= 0 or C0_Ce <= 0:
        raise ValueError(f"Initial concentrations must be positive: C0_MO={C0_MO}, C0_Ce={C0_Ce}")
    if MO_Dbulk <= 0 or CE_Dbulk <= 0:
        raise ValueError(f"Partition coefficients must be positive: MO_Dbulk={MO_Dbulk}, CE_Dbulk={CE_Dbulk}")
    
    # Check if U and Th parameters are available
    compute_u_th = (C0_U is not None and C0_Th is not None and 
                   U_Dbulk is not None and TH_Dbulk is not None and
                   C0_U > 0 and C0_Th > 0 and U_Dbulk > 0 and TH_Dbulk > 0)
    
    # Check if W parameters are available
    compute_w = (C0_W is not None and W_Dbulk is not None and
                 C0_W > 0 and W_Dbulk > 0)
    
    # Check if Nb and La parameters are available
    compute_nb_la = (C0_Nb is not None and C0_La is not None and
                     NB_Dbulk is not None and LA_Dbulk is not None and
                     C0_Nb > 0 and C0_La > 0 and NB_Dbulk > 0 and LA_Dbulk > 0)
    
    # Pre-compute R0 from delta0
    R0 = 1.0 + delta0 / 1000.0
    
    results = {}
    
    for alpha in alphas:
        rows = []
        validation_errors = []
        
        for F in F_list:
            row = {'F': F}
            
            if F == 0.0:
                # F=0: no fluid extracted yet
                # Residual rock = initial rock
                row['cs_Mo'] = C0_MO
                row['cs_Ce'] = C0_Ce
                row['fs'] = 1.0  # All Mo still in rock
                row['Rs'] = R0
                row['delta_s'] = delta0
                # Cumulative fluid at F=0 is undefined (no fluid extracted)
                row['cf_cum_Mo'] = np.nan
                row['cf_cum_Ce'] = np.nan
                row['MoCe_f_cum'] = np.nan
                row['Rf_cum'] = np.nan
                row['delta_f_cum'] = np.nan
                
                # U and Th at F=0
                if compute_u_th:
                    row['cs_U'] = C0_U
                    row['cs_Th'] = C0_Th
                    row['cf_cum_U'] = np.nan
                    row['cf_cum_Th'] = np.nan
                    row['UTh_f_cum'] = np.nan
                else:
                    row['cs_U'] = np.nan
                    row['cs_Th'] = np.nan
                    row['cf_cum_U'] = np.nan
                    row['cf_cum_Th'] = np.nan
                    row['UTh_f_cum'] = np.nan
                
                # W at F=0
                if compute_w:
                    row['cs_W'] = C0_W
                    row['cf_cum_W'] = np.nan
                    row['WTh_f_cum'] = np.nan
                    row['MoW_f_cum'] = np.nan
                    row['WCe_f_cum'] = np.nan
                else:
                    row['cs_W'] = np.nan
                    row['cf_cum_W'] = np.nan
                    row['WTh_f_cum'] = np.nan
                    row['MoW_f_cum'] = np.nan
                    row['WCe_f_cum'] = np.nan
                
                # Nb and La at F=0
                if compute_nb_la:
                    row['cs_Nb'] = C0_Nb
                    row['cs_La'] = C0_La
                    row['cf_cum_Nb'] = np.nan
                    row['cf_cum_La'] = np.nan
                    row['NbLa_f_cum'] = np.nan
                else:
                    row['cs_Nb'] = np.nan
                    row['cs_La'] = np.nan
                    row['cf_cum_Nb'] = np.nan
                    row['cf_cum_La'] = np.nan
                    row['NbLa_f_cum'] = np.nan
                
                # Ce/Mo at F=0 (inverse of Mo/Ce)
                row['CeMo_f_cum'] = np.nan
            else:
                # (1) Residual-rock concentrations (Chen Eq. 2)
                # cs_i(F) = c0_i * (1 - F)^((1/D_i) - 1)
                cs_Mo = C0_MO * ((1.0 - F) ** ((1.0 / MO_Dbulk) - 1.0))
                cs_Ce = C0_Ce * ((1.0 - F) ** ((1.0 / CE_Dbulk) - 1.0))
                
                row['cs_Mo'] = cs_Mo
                row['cs_Ce'] = cs_Ce
                
                # (2) Mo fraction remaining in rock (Chen Eq. 3)
                # fs(F) = (1 - F) * cs_Mo(F) / c0_Mo
                fs = (1.0 - F) * cs_Mo / C0_MO
                row['fs'] = fs
                
                # Validation: fs_simpl = (1 - F)^(1/D_Mo)
                fs_simpl = (1.0 - F) ** (1.0 / MO_Dbulk)
                if not np.isclose(fs, fs_simpl, rtol=1e-6):
                    validation_errors.append(f"F={F}: fs identity mismatch: fs={fs:.8f}, fs_simpl={fs_simpl:.8f}")
                
                # (3) Residual rock isotopes in ratio space
                # Rs(F) = R0 * fs(F)^(alpha - 1)
                Rs = R0 * (fs ** (alpha - 1.0))
                delta_s = (Rs - 1.0) * 1000.0
                
                row['Rs'] = Rs
                row['delta_s'] = delta_s
                
                # (4) CUMULATIVE extracted-fluid concentrations (mass balance)
                # cf_cum_i(F) = (c0_i - (1 - F)*cs_i(F)) / F
                cf_cum_Mo = (C0_MO - (1.0 - F) * cs_Mo) / F
                cf_cum_Ce = (C0_Ce - (1.0 - F) * cs_Ce) / F
                
                row['cf_cum_Mo'] = cf_cum_Mo
                row['cf_cum_Ce'] = cf_cum_Ce
                
                # Mass balance validation: c0_i ≈ (1-F)*cs_i + F*cf_cum_i
                mass_balance_Mo = (1.0 - F) * cs_Mo + F * cf_cum_Mo
                mass_balance_Ce = (1.0 - F) * cs_Ce + F * cf_cum_Ce
                if not np.isclose(C0_MO, mass_balance_Mo, rtol=1e-6):
                    validation_errors.append(f"F={F}: Mo mass balance error: C0={C0_MO}, reconstructed={mass_balance_Mo}")
                if not np.isclose(C0_Ce, mass_balance_Ce, rtol=1e-6):
                    validation_errors.append(f"F={F}: Ce mass balance error: C0={C0_Ce}, reconstructed={mass_balance_Ce}")
                
                # (5) CUMULATIVE extracted-fluid isotopes (closed form)
                # Rf_cum(F) = R0 * (1 - fs(F)^alpha) / (1 - fs(F))
                numerator = 1.0 - (fs ** alpha)
                denominator = 1.0 - fs
                
                if denominator == 0:
                    # This would only happen if fs=1, but fs<1 for F>0
                    Rf_cum = np.nan
                else:
                    Rf_cum = R0 * numerator / denominator
                
                delta_f_cum = (Rf_cum - 1.0) * 1000.0 if np.isfinite(Rf_cum) else np.nan
                
                row['Rf_cum'] = Rf_cum
                row['delta_f_cum'] = delta_f_cum
                
                # (6) Mo/Ce ratio in cumulative fluid
                if cf_cum_Ce != 0 and np.isfinite(cf_cum_Ce):
                    MoCe_f_cum = cf_cum_Mo / cf_cum_Ce
                else:
                    MoCe_f_cum = np.nan
                
                row['MoCe_f_cum'] = MoCe_f_cum
                
                # (7) U and Th calculations (same Rayleigh law)
                if compute_u_th:
                    # Residual rock concentrations for U and Th
                    cs_U = C0_U * ((1.0 - F) ** ((1.0 / U_Dbulk) - 1.0))
                    cs_Th = C0_Th * ((1.0 - F) ** ((1.0 / TH_Dbulk) - 1.0))
                    
                    row['cs_U'] = cs_U
                    row['cs_Th'] = cs_Th
                    
                    # Cumulative extracted-fluid concentrations for U and Th
                    cf_cum_U = (C0_U - (1.0 - F) * cs_U) / F
                    cf_cum_Th = (C0_Th - (1.0 - F) * cs_Th) / F
                    
                    row['cf_cum_U'] = cf_cum_U
                    row['cf_cum_Th'] = cf_cum_Th
                    
                    # Mass balance validation for U and Th
                    mass_balance_U = (1.0 - F) * cs_U + F * cf_cum_U
                    mass_balance_Th = (1.0 - F) * cs_Th + F * cf_cum_Th
                    if not np.isclose(C0_U, mass_balance_U, rtol=1e-6):
                        validation_errors.append(f"F={F}: U mass balance error: C0={C0_U}, reconstructed={mass_balance_U}")
                    if not np.isclose(C0_Th, mass_balance_Th, rtol=1e-6):
                        validation_errors.append(f"F={F}: Th mass balance error: C0={C0_Th}, reconstructed={mass_balance_Th}")
                    
                    # U/Th ratio in cumulative fluid
                    if cf_cum_Th != 0 and np.isfinite(cf_cum_Th) and np.isfinite(cf_cum_U):
                        UTh_f_cum = cf_cum_U / cf_cum_Th
                    else:
                        UTh_f_cum = np.nan
                    
                    row['UTh_f_cum'] = UTh_f_cum
                else:
                    row['cs_U'] = np.nan
                    row['cs_Th'] = np.nan
                    row['cf_cum_U'] = np.nan
                    row['cf_cum_Th'] = np.nan
                    row['UTh_f_cum'] = np.nan
                
                # (8) W calculations (same Rayleigh law)
                if compute_w:
                    # Residual rock concentration for W
                    cs_W = C0_W * ((1.0 - F) ** ((1.0 / W_Dbulk) - 1.0))
                    row['cs_W'] = cs_W
                    
                    # Cumulative extracted-fluid concentration for W
                    cf_cum_W = (C0_W - (1.0 - F) * cs_W) / F
                    row['cf_cum_W'] = cf_cum_W
                    
                    # Mass balance validation for W
                    mass_balance_W = (1.0 - F) * cs_W + F * cf_cum_W
                    if not np.isclose(C0_W, mass_balance_W, rtol=1e-6):
                        validation_errors.append(f"F={F}: W mass balance error: C0={C0_W}, reconstructed={mass_balance_W}")
                    
                    # W/Th ratio in cumulative fluid (requires both W and Th)
                    if compute_u_th and cf_cum_Th != 0 and np.isfinite(cf_cum_Th) and np.isfinite(cf_cum_W):
                        WTh_f_cum = cf_cum_W / cf_cum_Th
                    else:
                        WTh_f_cum = np.nan
                    row['WTh_f_cum'] = WTh_f_cum
                    
                    # Mo/W ratio in cumulative fluid
                    if cf_cum_W != 0 and np.isfinite(cf_cum_W) and np.isfinite(cf_cum_Mo):
                        MoW_f_cum = cf_cum_Mo / cf_cum_W
                    else:
                        MoW_f_cum = np.nan
                    row['MoW_f_cum'] = MoW_f_cum

                    # W/Ce ratio in cumulative fluid
                    if (cf_cum_Ce != 0 and np.isfinite(cf_cum_Ce) and
                            np.isfinite(cf_cum_W)):
                        WCe_f_cum = cf_cum_W / cf_cum_Ce
                    else:
                        WCe_f_cum = np.nan
                    row['WCe_f_cum'] = WCe_f_cum
                else:
                    row['cs_W'] = np.nan
                    row['cf_cum_W'] = np.nan
                    row['WTh_f_cum'] = np.nan
                    row['MoW_f_cum'] = np.nan
                    row['WCe_f_cum'] = np.nan
                
                # (9) Nb and La calculations (same Rayleigh law)
                if compute_nb_la:
                    cs_Nb = C0_Nb * ((1.0 - F) ** ((1.0 / NB_Dbulk) - 1.0))
                    cs_La = C0_La * ((1.0 - F) ** ((1.0 / LA_Dbulk) - 1.0))
                    row['cs_Nb'] = cs_Nb
                    row['cs_La'] = cs_La
                    
                    cf_cum_Nb = (C0_Nb - (1.0 - F) * cs_Nb) / F
                    cf_cum_La = (C0_La - (1.0 - F) * cs_La) / F
                    row['cf_cum_Nb'] = cf_cum_Nb
                    row['cf_cum_La'] = cf_cum_La
                    
                    if cf_cum_La != 0 and np.isfinite(cf_cum_La) and np.isfinite(cf_cum_Nb):
                        NbLa_f_cum = cf_cum_Nb / cf_cum_La
                    else:
                        NbLa_f_cum = np.nan
                    row['NbLa_f_cum'] = NbLa_f_cum
                else:
                    row['cs_Nb'] = np.nan
                    row['cs_La'] = np.nan
                    row['cf_cum_Nb'] = np.nan
                    row['cf_cum_La'] = np.nan
                    row['NbLa_f_cum'] = np.nan
                
                # (10) Ce/Mo ratio in cumulative fluid (inverse of Mo/Ce)
                if cf_cum_Mo != 0 and np.isfinite(cf_cum_Mo) and np.isfinite(cf_cum_Ce):
                    CeMo_f_cum = cf_cum_Ce / cf_cum_Mo
                else:
                    CeMo_f_cum = np.nan
                row['CeMo_f_cum'] = CeMo_f_cum
            
            rows.append(row)
        
        # Check for validation errors
        if validation_errors:
            # Return warning structure instead of raising
            print(f"Warning for alpha={alpha}: {len(validation_errors)} validation issues")
            for err in validation_errors[:5]:  # Show first 5
                print(f"  - {err}")
        
        # Create DataFrame
        df = pd.DataFrame(rows)
        
        # Reorder columns for clarity (include U, Th, and W)
        column_order = [
            'F', 'cs_Mo', 'cs_Ce', 'cs_U', 'cs_Th', 'cs_W', 'cs_Nb', 'cs_La',
            'fs', 'Rs', 'delta_s',
            'cf_cum_Mo', 'cf_cum_Ce', 'cf_cum_U', 'cf_cum_Th', 'cf_cum_W', 'cf_cum_Nb', 'cf_cum_La',
            'MoCe_f_cum', 'CeMo_f_cum', 'UTh_f_cum', 'WTh_f_cum', 'MoW_f_cum', 'WCe_f_cum', 'NbLa_f_cum',
            'Rf_cum', 'delta_f_cum'
        ]
        # Only include columns that exist
        column_order = [c for c in column_order if c in df.columns]
        df = df[column_order]
        
        results[alpha] = df

    return results


@st.cache_data(show_spinner=False)
def compute_filter_reequilibration(
    start: Dict[str, float],
    alpha_r_list: Tuple[float, ...] = ALPHA_R_DEFAULT,
    D_Mo_filter: float = None,
    D_W_filter: float = None,
    C_rockF_Mo: float = 3.96,
    C_rockF_W: float = 1.22,
    dlogfO2_end: float = None,
    iso_mode: str = "rayleigh",
) -> Dict[float, pd.DataFrame]:
    """
    Stage-2 REDOX FILTER — fluid–rock RE-EQUILIBRATION mechanism (D-controlled).

    Physical basis (reactive-transport / Ague-style filter)
    -------------------------------------------------------
    The proportional 'ratio' law (g = 10^(0.44·Δ)) assumes the fluid is
    Mo-SATURATED, but the Stage-1 fluid (0.5–2.7 ppm) is 600–1000× below the
    Bali solubility ceiling (~300–2700 ppm at slab conditions) — a ceiling that
    far above the load removes nothing. For dilute fluids the operative
    mechanism is PARTITIONING: the fluid exchanges Mo with the filter's mineral
    assemblage, buffered at the filter's ABSOLUTE fO₂ (ΔFMQ_filter, a
    petrological buffer floor: graphite/CCO-class metasediment ≈ FMQ −0.5…−4).
    At the filter's low fO₂ the bulk D is large, so a fully equilibrated fluid
    leaves with the rock-buffered value

        C_eq = C_rock^filter / D_bulk(T_seg, ΔFMQ_filter)     [rock-dominated]

    and the exchange progress f_ex (0 → 1) interpolates between no exchange and
    full equilibration of the MASS:

        C_f(f_ex)  = (1−f_ex)·C_start + f_ex·C_eq

    Isotopes (the fluid evolves heavier, never lighter):
    the Mo removed into the filter solids is the REDUCED species (Mo⁶⁺→Mo⁴⁺,
    weaker/longer bonds) and is therefore isotopically LIGHT relative to the
    fluid — the same polarity as Stage 1, where rutile (distorted octahedra)
    hosts light Mo and the tetrahedral aqueous molybdate is heavy (Chen 2019;
    Tang 2026). α_r ≡ R_fluid/R_precipitate > 1; magnitude justified by 1/T²
    scaling of measured low-T redox factors (0.8–3.0 ‰ at 25 °C → 0.06–0.27 ‰
    at 730–815 °C, i.e. α_r ≈ 1.0001–1.0003, coincident with the Stage-1
    family by scaling, not by assumption).

    Two isotope-evolution endmembers (`iso_mode`) share IDENTICAL elemental
    evolution (C_f, g, W, ratios, uptake mass) and differ ONLY in δ:

      'rayleigh' (open-system): each removed increment is isolated in the
          rock the moment it forms; integrating dδ = −1000·lnα_r·dln g gives
              δ_f(g) = δ_start − 1000·ln(α_r)·ln(g)     [unbounded as g→0]
          cumulative uptake δ = mass-balance complement
              δ_u = (δ_start − g·δ_f)/(1−g).
      'batch' (closed-system): the accumulated uptake stays in isotopic
          contact and shares one equilibrium with the fluid at the final g:
              δ_f(g) = δ_start + (1−g)·1000·ln(α_r)     [≤ +1000·lnα_r]
              δ_u    = δ_f − 1000·ln(α_r)
          (exact mass balance by construction: g·δ_f+(1−g)·δ_u = δ_start).

    Both endmembers make the surviving fluid HEAVIER; g = 1 ⇒ zero shift in
    both. The filter's PRE-EXISTING Mo inventory is assumed isotopically
    ISOLATED (refractory; no dissolution–reprecipitation on transit
    timescales) — full exchange with it would instead drag the fluid toward
    the sediment δ (≈ −1.8 ‰, lighter), a rejected endmember.
    Consistency: zero removal (clip, g = 1) ⇒ zero isotopic shift. W is
    treated analogously with its own filter D; fixed-D elements (Ce, Nb, Th,
    La) pass through unchanged. C_eq is clipped ≤ C_start (the filter never
    ADDS Mo). Removal % is EMERGENT: (1 − g)·100.

    Returns the standard Stage-2 schema (all plotting and tables consume it): 'g' = C_f/C_start (survival), removal_pct,
    dlogfO2 (constant, informational), is_endpoint (f_ex = 1), cf_cum_*,
    ratios, delta_f_cum, and the mass-balance ledger (Mo_precip_cum /
    delta_precip_cum = Mo and δ taken up by the filter rock).
    """
    if start is None:
        return {}
    Mo0 = start.get('Mo'); d0 = start.get('delta')
    if (Mo0 is None or d0 is None or not (Mo0 > 0)
            or not np.isfinite(Mo0) or not np.isfinite(d0)):
        return {}
    if D_Mo_filter is None or not np.isfinite(D_Mo_filter) or D_Mo_filter <= 0:
        return {}
    W0 = start.get('W'); Ce0 = start.get('Ce'); U0 = start.get('U')
    Th0 = start.get('Th'); Nb0 = start.get('Nb'); La0 = start.get('La')

    C_eq = min(C_rockF_Mo / D_Mo_filter, Mo0)          # clip: no Mo addition
    C_eqW = (min(C_rockF_W / D_W_filter, W0)
             if (W0 is not None and D_W_filter and D_W_filter > 0) else W0)

    def _ratio(num, den):
        if (num is None or den is None or not np.isfinite(num)
                or not np.isfinite(den) or den == 0):
            return np.nan
        return num / den

    NbLa_const = _ratio(Nb0, La0)
    # Progress points: regular f_ex grid PLUS rows at the exact marker survivals
    # (g = 0.75/0.50/0.25/0.10) where REACHABLE — so the removal-threshold
    # diamonds sit exactly on the curve. The curve is TRUNCATED at the
    # segment's equilibrium floor g_eq = C_eq/Mo0 (f_ex = 1, is_endpoint row):
    # ΔFMQ_filter sets that floor through D_bulk, so a deeper (more reducing)
    # filter draws a visibly longer curve. Survivals below g_eq are NOT drawn —
    # this filter cannot reach them.
    g_eq = C_eq / Mo0
    pts = [(round(0.05 * i, 2), None) for i in range(0, 21)]   # (f_ex, forced g)
    if g_eq < 1.0:
        for g_t in G_HIGHLIGHT:
            if g_eq <= g_t < 1.0:
                pts.append(((1.0 - g_t) / (1.0 - g_eq), g_t))
    pts = sorted(set(pts), key=lambda x: x[0])

    _is_batch = (str(iso_mode).strip().lower() == "batch")
    results: Dict[float, pd.DataFrame] = {}
    for alpha_r in alpha_r_list:
        ln_ar = np.log(alpha_r)
        rows = []
        for fex, g_forced in pts:
            cf = (Mo0 * g_forced) if g_forced is not None \
                else (1.0 - fex) * Mo0 + fex * C_eq
            cfW = ((1.0 - fex) * W0 + fex * C_eqW) if W0 is not None else None
            g = g_forced if g_forced is not None else cf / Mo0
            # Removed increments are LIGHT (α_r = R_fluid/R_precip > 1) →
            # fluid drifts HEAVIER in BOTH endmembers; g=1 → zero shift.
            if g >= 1.0:
                dl = d0
                d_rock_uptake = d0 - 1000.0 * ln_ar   # first-increment limit
            elif g > 0:
                if _is_batch:
                    # closed-system: one equilibrium at the final g
                    dl = d0 + (1.0 - g) * 1000.0 * ln_ar
                    d_rock_uptake = dl - 1000.0 * ln_ar
                else:
                    # open-system Rayleigh on the survival fraction g
                    dl = d0 - 1000.0 * ln_ar * np.log(g)
                    d_rock_uptake = ((Mo0 * d0) - (cf * dl)) / (Mo0 - cf) \
                        if (Mo0 - cf) > 1e-12 else np.nan
            else:
                dl = np.nan
                d_rock_uptake = np.nan
            removed = Mo0 - cf
            row = {
                'g': g, 'fex': fex, 'alpha_r': alpha_r,
                'isotope_mode': 'batch' if _is_batch else 'rayleigh',
                'dlogfO2': dlogfO2_end if dlogfO2_end is not None else np.nan,
                'removal_pct': (1.0 - g) * 100.0,
                'g_W': (cfW / W0) if (W0 not in (None, 0) and cfW is not None) else np.nan,
                'g_U': 1.0,
                'is_endpoint': int(np.isclose(fex, 1.0)),
                'Mo_start': Mo0, 'Mo_eq_unclipped': C_rockF_Mo / D_Mo_filter,
                'Mo_eq': C_eq, 'D_Mo_filter': D_Mo_filter,
                'C_rockF_Mo': C_rockF_Mo, 'delta_start': d0,
                'W_start': W0 if W0 is not None else np.nan,
                'W_eq': C_eqW if C_eqW is not None else np.nan,
                'W_eq_unclipped': (
                    C_rockF_W / D_W_filter
                    if D_W_filter is not None and np.isfinite(D_W_filter) and D_W_filter > 0
                    else np.nan
                ),
                'D_W_filter': D_W_filter if D_W_filter is not None else np.nan,
                'C_rockF_W': C_rockF_W,
                'cf_cum_Mo': cf, 'cf_cum_Ce': Ce0, 'cf_cum_U': U0,
                'cf_cum_Th': Th0, 'cf_cum_W': cfW, 'cf_cum_Nb': Nb0,
                'cf_cum_La': La0,
                'MoCe_f_cum': _ratio(cf, Ce0), 'CeMo_f_cum': _ratio(Ce0, cf),
                'MoW_f_cum': _ratio(cf, cfW), 'UTh_f_cum': _ratio(U0, Th0),
                'WTh_f_cum': _ratio(cfW, Th0), 'NbLa_f_cum': NbLa_const,
                'WCe_f_cum': _ratio(cfW, Ce0),
                'Rf_cum': (1.0 + dl / 1000.0) if np.isfinite(dl) else np.nan,
                'delta_f_cum': dl,
                'Mo_precip_cum': removed,
                'delta_precip_cum': d_rock_uptake,
            }
            rows.append(row)
        results[alpha_r] = pd.DataFrame(rows)
    return results


def extract_reduction_start(
    results_stage1: Dict[float, pd.DataFrame],
    F_branch: float = F_BRANCH,
    alpha_ref: float = ALPHA_EXT_REF,
) -> Optional[Dict[str, float]]:
    """
    Read the Stage-1 cumulative extracted-fluid composition at F = F_branch.

    This is the START point for the Stage-2 transport-reduction Rayleigh. The Mo
    concentration and δ⁹⁸/⁹⁵Mo are taken from the extraction curve at F_branch using
    `alpha_ref` — a fixed-reference convention (α = 1.0003) giving all Stage-2
    trajectories a common origin. The suppressed Stage-1 α-spread of δ_start is
    negligible when Mo is efficiently extracted but grows to ±0.1 ‰ for
    Mo-retentive segments (f_s → 1); see Supplementary §6.5. If alpha_ref is not
    present, the closest available α is used. The non-Mo elements (Ce, U, Th, W,
    Nb, La) are held constant through Stage 2.

    Returns
    -------
    dict with keys Mo, delta, Ce, U, Th, W, Nb, La (floats or None), or None if
    the branch point cannot be located.
    """
    if not results_stage1:
        return None
    # Pick the reference α if available, else the α closest to alpha_ref
    if alpha_ref in results_stage1:
        df = results_stage1[alpha_ref]
    else:
        _keys = sorted(results_stage1.keys())
        df = results_stage1[min(_keys, key=lambda a: abs(a - alpha_ref))]
    if df is None or df.empty or 'F' not in df.columns:
        return None
    match = df[np.isclose(df['F'].astype(float), F_branch, atol=1e-9)]
    if match.empty:
        return None
    r = match.iloc[0]

    def _g(col):
        if col in df.columns and pd.notna(r[col]):
            return float(r[col])
        return None

    Mo_start = _g('cf_cum_Mo')
    delta_start = _g('delta_f_cum')
    if Mo_start is None or delta_start is None or Mo_start <= 0:
        return None
    return {
        'Mo': Mo_start,
        'delta': delta_start,
        'Ce': _g('cf_cum_Ce'),
        'U': _g('cf_cum_U'),
        'Th': _g('cf_cum_Th'),
        'W': _g('cf_cum_W'),
        'Nb': _g('cf_cum_Nb'),
        'La': _g('cf_cum_La'),
    }


# ============================================================================
# UI RENDERING FUNCTIONS
# ============================================================================

def create_fractionation_plot_vs_mo(
    title: str = "δ⁹⁸/⁹⁵Mo vs Mo (Cumulative Extracted Fluid)",
    log_x: bool = True,
    arc_df: Optional[pd.DataFrame] = None,
    style_df: Optional[pd.DataFrame] = None,
    arc_settings: Optional[Dict] = None,
    dfmq: float = None,
    mode_rut_pct: float = None
) -> go.Figure:
    """Create plot of δ⁹⁸/⁹⁵Mo vs Mo concentration (cumulative fluid) for all alphas.

    Plot H: X = Mo (µg/g) in the cumulative extracted fluid (the variable, not a
    ratio), Y = δ⁹⁸/⁹⁵Mo in the cumulative fluid. Mirrors Plot G in style.
    """

    fig = go.Figure()


    if arc_df is not None and not arc_df.empty:
        if arc_settings is None:
            arc_settings = {'color_by': 'Arc', 'symbol_by': 'Arc', 'show_error_bars': True,
                            'show_labels': False, 'legend_by': 'Arc', 'filters': {}}
        fig = add_arc_overlay_to_figure(
            fig=fig, arc_df=arc_df, style_df=style_df,
            x_col='Mo_ppm', y_col='d98Mo',
            x_err_col=None,
            y_err_col='d98Mo_2SE' if arc_settings.get('show_error_bars', True) else None,
            color_by=arc_settings.get('color_by', 'Arc'),
            symbol_by=arc_settings.get('symbol_by', 'Arc'),
            show_error_bars=arc_settings.get('show_error_bars', True),
            show_labels=arc_settings.get('show_labels', False),
            legend_by=arc_settings.get('legend_by', 'Arc'),
            filters=arc_settings.get('filters', {})
        )


    fig.update_layout(
        template='simple_white', width=700, height=480, margin=dict(l=60, r=250, t=30, b=55),
        paper_bgcolor='white', plot_bgcolor='white', font=dict(color='#222222'),
        showlegend=False, hovermode='closest',
        xaxis=dict(title='Mo (µg/g, cumulative fluid)', type='log' if log_x else 'linear', autorange=True,
                   showgrid=False, showline=True, mirror=True, linecolor='black', linewidth=1,
                   ticks='outside', dtick=1 if log_x else None, tickformat='.3~g' if log_x else None,
                   minor=dict(ticks='outside', showgrid=False)),
        yaxis=dict(title='δ⁹⁸/⁹⁵Mo (‰, cumulative fluid)', autorange=True, showgrid=False, showline=True,
                   mirror=True, linecolor='black', linewidth=1, ticks='outside', zeroline=False, dtick=0.2,
                   minor=dict(ticks='outside', showgrid=False)),
    )
    if title:
        # conditions block to the right of the plot box
        fig.add_annotation(xref='paper', yref='paper', x=1.03, y=0.98, xanchor='left', yanchor='top',
                           text=title, showarrow=False, align='left', font=dict(size=11, color='#222222'))

    return fig


# ── Stage colour scheme ─────────────────────────────────────────────────────
# Extraction (Stage 1, fluid GAINS Mo) = RED ; Reduction (Stage 2, fluid LOSES
# Mo) = BLUE. Line style encodes protolith: Unmodified = solid, Modified = dashed.
_STAGE_COLOR_EXTRACT = '#DC143C'   # red   (crimson)
_STAGE_COLOR_REDUCE = '#1F5FBF'    # blue

# Fallback colour cycle for segments outside the named regions
# (Nicaragua = '#0072B2', Costa Rica = '#D55E00') — Okabe–Ito, minus those two.
_SEGMENT_COLORS = ['#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#000000']


def build_plot_line_table(fig: go.Figure) -> pd.DataFrame:
    """Return paired numeric X/Y columns from the visible model traces.

    Preserve trace order and point order, including separate coordinates for
    every alpha. Missing trailing points are blank, never interpolated.
    """
    columns = {}
    for trace in fig.data:
        meta = trace.meta
        if not isinstance(meta, dict) or not meta.get('part2_model_line'):
            continue
        if trace.visible is False or trace.visible == 'legendonly':
            continue
        curve = len(columns) // 3 + 1
        suffix = f"curve{curve}_alpha{meta['alpha']}"
        columns[f"X_{meta['x_column']}_{suffix}"] = pd.Series(trace.x, dtype=float)
        columns[f"Y_delta98Mo_permil_{suffix}"] = pd.Series(trace.y, dtype=float)
        columns[f"{meta['progress_column']}_{suffix}"] = pd.Series(
            np.asarray(trace.customdata)[:, 0], dtype=float)
    return pd.DataFrame(columns)


def _add_stage_curves(fig, results, x_col, color, dash, label_prefix,
                      progress_col, marker_vals, y_col='delta_f_cum', opacity=1.0):
    """Draw one Rayleigh-stage family (extraction OR reduction) onto `fig`.

    One trace per α (each carries its α in the legend name AND the hover), sharing
    `color` + `dash`. `progress_col` is 'F' (extraction) or 'g' (reduction);
    `marker_vals` are the progress values that get a labelled marker
    (F ∈ {0.1,0.2,0.3} for extraction, g = 0.25 for reduction).
    """
    if not results:
        return fig
    is_reduction = (progress_col == 'g')
    for a in sorted(results.keys()):
        df = results[a]
        if df is None or progress_col not in df.columns or x_col not in df.columns:
            continue
        dfp = df[df[x_col].notna() & (df[x_col] > 0) & df[y_col].notna()].copy()
        if dfp.empty:
            continue
        name = f"{label_prefix} · α={a}"
        _plabel = 'g' if is_reduction else 'F'
        fig.add_trace(go.Scatter(
            x=dfp[x_col], y=dfp[y_col], mode='lines',
            name=name, legendgroup=name,
            line=dict(color=color, width=1.1, dash=dash), opacity=opacity,
            hovertemplate=(
                f'{name}<br>{_plabel} = %{{customdata[0]:.2f}}<br>'
                'x = %{x:.4f}<br>δ⁹⁸/⁹⁵Mo = %{y:.3f}‰<extra></extra>'
            ),
            customdata=dfp[[progress_col]].values,
            meta={'part2_model_line': True, 'alpha': a, 'x_column': x_col,
                  'progress_column': progress_col},
            showlegend=True,
        ))
        hl = dfp[dfp[progress_col].isin(marker_vals)]
        if not hl.empty:
            if is_reduction:
                # Removal thresholds 25/50/75/90 % (g = 0.75/0.50/0.25/0.10).
                # Only 75 % carries the full "Mo lost" label; others just the %.
                mtext = [('75% Mo lost' if np.isclose(v, 0.25)
                          else f'{int(round((1.0 - v) * 100))}%')
                         for v in hl[progress_col]]
                msym = 'diamond'
                msize = 7
            else:
                mtext = [f'{int(round(v * 100))}%' for v in hl[progress_col]]
                msym = 'circle'
                msize = 6
            fig.add_trace(go.Scatter(
                x=hl[x_col], y=hl[y_col], mode='markers+text',
                name=f'{name} marker', legendgroup=name, showlegend=False,
                marker=dict(color=color, size=msize, symbol=msym,
                            line=dict(color='black', width=0.8)),
                text=mtext, textposition='top center',
                textfont=dict(size=9, color=color, family='Arial'),
                customdata=hl[[progress_col]].values, opacity=opacity,
                hovertemplate=(
                    f'{name}<br>{_plabel} = %{{customdata[0]:.2f}}<br>'
                    'x = %{x:.4f}<br>δ⁹⁸/⁹⁵Mo = %{y:.3f}‰<extra></extra>'
                ),
            ))
        # Reduction progress points are uniform diamonds; the curve itself is
        # truncated at the filter's equilibrium floor (is_endpoint row), so the
        # reachable extent — set by ΔFMQ_filter — is what the curve length shows.
    return fig


@st.fragment
def render_fluid_fractionation_section(
    C0_MO: float = None,
    C0_Ce: float = None,
    MO_Dbulk: float = None,
    CE_Dbulk: float = None,
    salinity: float = None,
    delta0: float = -0.21,
    alphas: Tuple[float, ...] = (1.0002, 1.0003, 1.0004),
    arc_df: Optional[pd.DataFrame] = None,
    style_df: Optional[pd.DataFrame] = None,
    arc_settings: Optional[Dict] = None
) -> None:
    """
    Render the complete fluid fractionation UI section in Streamlit.
    
    If parameters are not provided, attempts to read from st.session_state.
    
    Parameters
    ----------
    C0_MO, C0_Ce : float, optional
        Initial concentrations (ppm). If None, reads from session state.
    MO_Dbulk, CE_Dbulk : float, optional
        Bulk partition coefficients. If None, reads from session state.
    salinity : float, optional
        Selected salinity for display purposes.
    delta0 : float
        Initial δ⁹⁸/⁹⁵Mo value (‰)
    alphas : tuple of float
        Fractionation factors
    arc_df : pd.DataFrame, optional
        Arc sample data for overlay
    style_df : pd.DataFrame, optional
        Style definitions for arc data
    arc_settings : dict, optional
        Arc display settings
    
    Notes
    -----
    **Salinity Selection:**
    The bulk partition coefficients (D_bulk) for Mo, Ce, U, and Th are salinity-dependent.
    - Salinity is expressed in weight percent (wt%) NaCl
    - Internally, the Bali model converts wt% to molality (mol NaCl / kg H₂O) for calculations
    - Higher salinity → different D_bulk values → different fractionation curves
    
    **Source Parameters:**
    - C₀(Mo), C₀(Ce), C₀(U), C₀(Th): Initial concentrations from the selected source rock
    - D_bulk values: Bulk partition coefficients that depend on the selected salinity
    """
    pass  # section title is given by the dashboard expander

    # Driven entirely by the 'Part2_Segments' sheet.
    p2_segments = st.session_state.get('p2_segments') or []
    if not p2_segments:
        st.warning(
            "**No slab segments found.**\n\n"
            "Add a sheet named **'Part2_Segments'** to the input workbook (one row per slab "
            "segment: segment, P_bar, T_C, T_K, phi, ΔFMQ (log units), mode_grt/cpx/rut, "
            "C0_MO … C0_La, d98Mo, d98Mo 2SD, rock_id)."
        )
        return

    # Available salinities come from a ready segment's own D_bulk table (engine fixed grid).
    _ready_seg = next((s for s in p2_segments
                       if s.get('ready') and s.get('results_df') is not None), None)
    if _ready_seg is not None:
        available_salinities = sorted(_ready_seg['results_df']['NaCl_wt_pct'].unique().tolist())
    else:
        available_salinities = [0.001, 5.0, 10.0, 15.0, 20.0]
    
    # Salinity selection with radio-button style behavior
    st.markdown("#### Fluid salinity (wt.% NaCl)")
    
    # Initialize session state for fluid fractionation salinity - default to 5%
    if 'ff_selected_salinity' not in st.session_state:
        st.session_state.ff_selected_salinity = 5.0 if 5.0 in available_salinities else available_salinities[0]
    
    # Salinity labels for display
    salinity_labels = {
        0.001: "0.001%",
        5.0: "5%",
        10.0: "10%",
        15.0: "15%",
        20.0: "20%"
    }
    
    # Use radio buttons for single-click selection/deselection
    salinity_options = [salinity_labels.get(sal, f"{sal}%") for sal in available_salinities]
    salinity_map = {salinity_labels.get(sal, f"{sal}%"): sal for sal in available_salinities}
    
    current_label = salinity_labels.get(st.session_state.ff_selected_salinity, f"{st.session_state.ff_selected_salinity}%")
    current_index = salinity_options.index(current_label) if current_label in salinity_options else 0
    
    selected_label = st.radio(
        "Select salinity:",
        options=salinity_options,
        index=current_index,
        horizontal=True,
        key="ff_salinity_radio",
        label_visibility="collapsed"
    )
    
    new_salinity = salinity_map[selected_label]
    st.session_state.ff_selected_salinity = new_salinity
    selected_salinity = new_salinity
    
    # Toggle the φ (Zr-in-rutile) slab correction.
    with st.expander("Model settings — φ correction, curve visibility, Stage-2 filter and isotope law", expanded=False):
        # φ is read per segment from the 'Part2_Segments' sheet (column `phi`).
        apply_phi600 = st.checkbox(
            "Apply φ (Zr-in-rutile) slab correction — per-segment φ from Excel (Chen et al. 2019)",
            value=True,
            key="ff_apply_phi600",
            help=(
                "ON: D = D(bracket at the segment's T_K) × φ(segment) — the φ value is read "
                "from the 'Part2_Segments' sheet (column phi), computed from each segment's T_K. "
                "OFF: uses the bracket coefficients without the φ mineral correction. "
                )
        )

        # ── Two-stage / two-protolith display toggles ──────────────────────────
        # Colour scheme (fixed): Extraction = red, Reduction = blue.
        # Line style: Unmodified protolith = solid, Modified protolith = dashed.
        st.markdown("#### Curve visibility")
        _vis_c1, _vis_c2 = st.columns(2)
        with _vis_c1:
            show_unmodified = st.checkbox(
                "Show **Unmodified** protolith (solid)",
                value=True, key="p2_show_unmodified",
                help="Pure-source Rayleigh paths (no sediment), drawn solid."
            )
            show_modified = st.checkbox(
                "Show **Modified** protolith (dashed)",
                value=True, key="p2_show_modified",
                help="Sediment-modified Rayleigh paths (uses the sediment sliders), drawn dashed."
            )
        with _vis_c2:
            show_extraction = st.checkbox(
                "Show **Extraction** stage (red)",
                value=True, key="p2_show_extraction",
                help="Stage 1 — slab-dehydration Rayleigh; the fluid GAINS Mo."
            )
            show_reduction = st.checkbox(
                "Show **Reduction** stage (blue)",
                value=True, key="p2_show_reduction",
                help="Stage 2 — redox filter; the fluid LOSES Mo by re-equilibrating "
                     "with the filter rock at its buffer fO₂ (C_eq = C_rockF/D(T, "
                     "ΔFMQ_filter), clipped ≤ start — removal EMERGENT, driven by "
                     "ΔFMQ_filter in the Part2_Segments sheet). Branches off "
                     f"the extraction curve at F = {F_BRANCH:.0%} "
                     f"(extraction α = {ALPHA_EXT_REF})."
            )

        # ── Filter mechanism: fluid–rock RE-EQUILIBRATION (the only Stage-2 law) ──
        # The fluid exchanges Mo with the filter rock at the filter's ABSOLUTE
        # buffer fO₂ (ΔFMQ_filter, e.g. graphite/CCO-class ≈ −0.5…−4); exit value is
        # rock-buffered, C_eq = C_rock/D_bulk(T, ΔFMQ_filter), clipped ≤ C_start (the
        # sediment never ADDS net Mo). The earlier proportional/saturation law was
        # (the fluid is ×600–1000 undersaturated, so a saturation law does not apply).
        rockF_Mo = st.number_input(
            "Filter rock Mo (ppm) — fallback for blank C0_MO_filter cells", min_value=0.0, max_value=50.0,
            value=3.96, step=0.01, format="%.2f", key="p2_rockF_Mo",
            help="Mo content of the filter's sediment (default: hemipelagic "
                 "endmember 3.96 ppm — the same measured endmember used for the "
                 "protolith mixing; treat as data, not a fit parameter). "
                 "Per-segment values in the 'C0_MO_filter' column of the "
                 "Part2_Segments sheet take precedence; this input is used only "
                 "where that cell is blank. W is "
                 "scaled proportionally from the hemipelagic W = 1.22 ppm. "
                 "Note: the filter rock's δ⁹⁸/⁹⁵Mo does NOT enter the model — its "
                 "pre-existing Mo is isotopically isolated (refractory); only the "
                 "newly removed increments fractionate (light into the solid, "
                 "fluid drifts heavier; α_r = R_fluid/R_precipitate).")
        rockF_W = 1.22 * (rockF_Mo / 3.96) if rockF_Mo > 0 else 1.22

        # ── Isotope evolution endmember (elemental evolution IDENTICAL in both) ──
        _iso_label = st.radio(
            "Stage-2 isotope evolution",
            ["Rayleigh (open-system)", "Batch (closed-system)"],
            index=1, horizontal=True, key="p2_iso_mode",
            help="Both endmembers share the SAME mass law (re-equilibration floor) "
                 "and identical Mo, g, Mo/W, Mo/Ce and uptake mass — they differ "
                 "only in δ⁹⁸/⁹⁵Mo. Rayleigh: removed increments isolated as they form "
                 "→ δ_f = δ_start − 1000·lnα_r·ln g (unbounded). Batch: accumulated "
                 "uptake stays in isotopic contact, one equilibrium at exit → "
                 "δ_f = δ_start + (1−g)·1000·lnα_r (capped at +1000·lnα_r). The "
                 "surviving fluid drifts HEAVIER in both (removed reduced Mo is "
                 "light; α_r = R_fluid/R_precipitate > 1)."
        )
        iso_mode = "batch" if _iso_label.startswith("Batch") else "rayleigh"

        # Custom α_r (exploration): 0 = off; any value >1 adds an extra curve.
        custom_ar = st.number_input(
            "Custom α_r (0 = off; e.g. 1.00008 for near-zero fractionation)",
            min_value=0.0, max_value=1.01, value=0.0, step=0.00001, format="%.5f",
            key="p2_custom_ar",
            help="Adds one extra reduction curve with this fractionation factor for "
                 "both protoliths (labelled 'custom'). The same-slab/Nicaragua-filter "
                 "scenario requires α_r ≈ 1.0001 under Rayleigh removal."
        )
        custom_ar = float(custom_ar) if custom_ar and custom_ar > 1.0 else None
        alpha_r_active_list = tuple(sorted(set(ALPHA_R_DEFAULT)
                                           | ({custom_ar} if custom_ar else set())))

    # ── Lava data overlay ──
    # Same lava data and colours as the Mo-W Thermodynamic model panels (data/arc_data.xlsx).
    arc_df, _style_file, _arc_errs = load_arc_data_from_file()
    _ARC_COL = {"Nicaragua": ("#fca9ff", "#9c329a"), "CostaRica": ("#93e65d", "#3a7a1e")}
    _rows = []
    if arc_df is not None and not arc_df.empty:
        _grp_opacity = {}
        if _style_file is not None and not _style_file.empty:
            for _, _r in _style_file[_style_file["KeyType"].astype(str) == "Group"].iterrows():
                _grp_opacity[str(_r["Key"])] = float(_r.get("Opacity", 1.0))
        for (_arc, _grp), _ in arc_df.groupby(["Arc", "Group"], dropna=False):
            _fill, _line = _ARC_COL.get(str(_arc), ("#a0a0a0", "#555555"))
            _rows.append({"KeyType": "Group", "Key": _grp, "Color": _fill, "Symbol": "circle", "Size": 9,
                          "Opacity": _grp_opacity.get(str(_grp), 1.0), "OutlineColor": _line,
                          "OutlineWidth": 1, "Visible": 1})
    style_df = pd.DataFrame(_rows) if _rows else None
    arc_settings = {"color_by": "Group", "symbol_by": "Group", "show_error_bars": True,
                    "show_labels": False, "legend_by": "Arc", "filters": {}}

    # ── Sediment-modified protolith sliders ──
    _p2s_col1, _p2s_col2 = st.columns(2)
    with _p2s_col1:
        sed_pct_val = st.slider(
            "Sediment fraction in slab (%)", 0, 100,
            value=5, step=1,
            key="p2_sed_pct_slider"
        )
    with _p2s_col2:
        carb_hemi_pct_val = st.slider(
            "Carb ← Sed type → Hemi (%)", 0, 100,
            value=50, step=1,
            key="p2_carb_hemi_slider",
            help="0 = pure Carbonates, 100 = pure Hemipelagic"
        )

    # ====================================================================
    # SOURCE SELECTION (single active row) + Fractionation factor α
    # ====================================================================
    salinity = selected_salinity

    def _seg_region(nm):
        n = str(nm).lower()
        return "Nicaragua" if "nicaragua" in n else ("Costa Rica" if "costa" in n else str(nm))

    # ── Single active slab segment (checkbox table) ──
    _p2src_rows = []
    for _s in p2_segments:
        _c0 = _s.get('C0', {}) or {}
        _dfmq_f = _s.get('dFMQ_filter')
        _dfmq_s = _s.get('dFMQ')
        _p2src_rows.append({
            "Select": False,
            "row": _s.get('row'), "segment": _s['name'], "protolith": _s.get('protolith'),
            "T_K": _s.get('T_K'), "phi": (round(_s['phi'], 4) if _s.get('phi') is not None else None),
            "ΔFMQ_slab": _dfmq_s,
            "ΔFMQ_filter": (_dfmq_f if (_dfmq_f is not None and np.isfinite(_dfmq_f)) else None),
            "Δlog fO₂": (round(_dfmq_f - _dfmq_s, 3)
                         if (_dfmq_f is not None and _dfmq_s is not None
                             and np.isfinite(_dfmq_f) and np.isfinite(_dfmq_s)) else None),
            "P_bar": _s.get('P_bar'), "mode_rut%": _s.get('mode_rut'),
            "C0_MO": _c0.get('MO'), "C0_Ce": _c0.get('Ce'), "C0_W": _c0.get('W'),
            "C0_U": _c0.get('U'), "C0_Th": _c0.get('Th'), "C0_Nb": _c0.get('Nb'), "C0_La": _c0.get('La'),
            "d98Mo": _s.get('d98Mo'), "d98Mo_2SD": _s.get('d98Mo_2SD'), "ready": _s.get('ready'),
        })
    _p2_order = [r["row"] for r in _p2src_rows]
    st.session_state["_p2_src_order"] = _p2_order
    if st.session_state.get("p2_active_row") not in _p2_order:
        # default active = first READY row, else first row
        _first_ready = next((r["row"] for r in _p2src_rows if r["ready"]), None)
        st.session_state["p2_active_row"] = _first_ready if _first_ready is not None else (_p2_order[0] if _p2_order else None)
    for _r in _p2src_rows:
        _r["Select"] = (_r["row"] == st.session_state.get("p2_active_row"))

    def _p2_src_on_select():
        _ed = st.session_state.get("p2_src_editor", {})
        _er = _ed.get("edited_rows", {}) if isinstance(_ed, dict) else {}
        _order = st.session_state.get("_p2_src_order", [])
        _newly = [int(p) for p, ch in _er.items() if isinstance(ch, dict) and ch.get("Select") is True]
        if _newly and _order and _newly[-1] < len(_order):
            st.session_state["p2_active_row"] = _order[_newly[-1]]
        if isinstance(_ed, dict):
            _ed["edited_rows"] = {}

    st.markdown("**Slab segment — sheet `Part2_Segments` of `data/sample_input_data.xlsx`**")
    st.caption("Tick one row to select the segment used for the calculation. "
               "Edit or add rows in the workbook to model other segments.")
    st.data_editor(
        pd.DataFrame(_p2src_rows),
        key="p2_src_editor",
        on_change=_p2_src_on_select,
        column_config={"Select": st.column_config.CheckboxColumn(
            "✓", help="Active slab segment", default=False, width="small")},
        disabled=[c for c in _p2src_rows[0].keys() if c != "Select"] if _p2src_rows else [],
        hide_index=True,
        use_container_width=True,
        height=180,
    )
    _visible = {st.session_state.get("p2_active_row")}

    # ── Fractionation factor α — Extraction (Stage 1) vs Reduction (Stage 2) ──
    # Defaults: only extraction α=1.0003 and reduction α_r=1.0003 are ON, per
    # protolith (→ 4 curves visible by default). The other α are computed for the
    # tables but hidden until toggled here.
    _alpha_col, _plot_col = st.columns([1, 1])
    _alpha_ctx = _alpha_col.container()
    _alpha_ctx.__enter__()
    st.markdown("**Extraction α (Stage 1) — red curves**")
    _ctrl_a, _ctrl_b = st.columns(2)
    with _ctrl_a:
        st.markdown("*Unmodified source (solid)*")
        _sa1, _sa2, _sa3 = st.columns([1, 1, 2])
        with _sa1:
            if st.button("All", key="alpha_src_all_btn", use_container_width=True):
                for _k in ('alpha_src_0002', 'alpha_src_0003', 'alpha_src_0004'):
                    st.session_state[_k] = True
        with _sa2:
            if st.button("None", key="alpha_src_none_btn", use_container_width=True):
                for _k in ('alpha_src_0002', 'alpha_src_0003', 'alpha_src_0004'):
                    st.session_state[_k] = False
        show_alpha_0002 = st.checkbox("α = 1.0002", value=False, key="alpha_src_0002")
        show_alpha_0003 = st.checkbox("α = 1.0003", value=True, key="alpha_src_0003")
        show_alpha_0004 = st.checkbox("α = 1.0004", value=False, key="alpha_src_0004")
    with _ctrl_b:
        st.markdown("*Sediment-modified source (dashed)*")
        _sd1, _sd2, _sd3 = st.columns([1, 1, 2])
        with _sd1:
            if st.button("All", key="alpha_sed_all_btn", use_container_width=True):
                for _k in ('alpha_sed_0002', 'alpha_sed_0003', 'alpha_sed_0004'):
                    st.session_state[_k] = True
        with _sd2:
            if st.button("None", key="alpha_sed_none_btn", use_container_width=True):
                for _k in ('alpha_sed_0002', 'alpha_sed_0003', 'alpha_sed_0004'):
                    st.session_state[_k] = False
        show_alpha_sed_0002 = st.checkbox("α = 1.0002", value=False, key="alpha_sed_0002")
        show_alpha_sed_0003 = st.checkbox("α = 1.0003", value=True, key="alpha_sed_0003")
        show_alpha_sed_0004 = st.checkbox("α = 1.0004", value=False, key="alpha_sed_0004")

    st.markdown("**Reduction α_r (Stage 2, redox filter) — blue curves**")
    _ctrl_c, _ctrl_d = st.columns(2)
    _AR_SRC_KEYS = ('alpha_red_src_0001', 'alpha_red_src_0002',
                    'alpha_red_src_0003', 'alpha_red_src_0004')
    _AR_SED_KEYS = ('alpha_red_sed_0001', 'alpha_red_sed_0002',
                    'alpha_red_sed_0003', 'alpha_red_sed_0004')
    with _ctrl_c:
        st.markdown("*Unmodified source (solid)*")
        _ra1, _ra2, _ra3 = st.columns([1, 1, 2])
        with _ra1:
            if st.button("All", key="alpha_red_src_all_btn", use_container_width=True):
                for _k in _AR_SRC_KEYS:
                    st.session_state[_k] = True
        with _ra2:
            if st.button("None", key="alpha_red_src_none_btn", use_container_width=True):
                for _k in _AR_SRC_KEYS:
                    st.session_state[_k] = False
        show_ar_src_0001 = st.checkbox("α_r = 1.00002 (near-zero)", value=False, key="alpha_red_src_0001")
        show_ar_src_0002 = st.checkbox("α_r = 1.0002", value=False, key="alpha_red_src_0002")
        show_ar_src_0003 = st.checkbox("α_r = 1.0003", value=True, key="alpha_red_src_0003")
        show_ar_src_0004 = st.checkbox("α_r = 1.0004", value=False, key="alpha_red_src_0004")
    with _ctrl_d:
        st.markdown("*Sediment-modified source (dashed)*")
        _rd1, _rd2, _rd3 = st.columns([1, 1, 2])
        with _rd1:
            if st.button("All", key="alpha_red_sed_all_btn", use_container_width=True):
                for _k in _AR_SED_KEYS:
                    st.session_state[_k] = True
        with _rd2:
            if st.button("None", key="alpha_red_sed_none_btn", use_container_width=True):
                for _k in _AR_SED_KEYS:
                    st.session_state[_k] = False
        show_ar_sed_0001 = st.checkbox("α_r = 1.00002 (near-zero)", value=False, key="alpha_red_sed_0001")
        show_ar_sed_0002 = st.checkbox("α_r = 1.0002", value=False, key="alpha_red_sed_0002")
        show_ar_sed_0003 = st.checkbox("α_r = 1.0003", value=True, key="alpha_red_sed_0003")
        show_ar_sed_0004 = st.checkbox("α_r = 1.0004", value=False, key="alpha_red_sed_0004")
    _alpha_ctx.__exit__(None, None, None)

    # ── α active sets (extraction: unmodified 'src' / modified 'sed') ──
    _active_alphas_src = {a for a, on in [(1.0002, show_alpha_0002), (1.0003, show_alpha_0003),
                                          (1.0004, show_alpha_0004)] if on}
    _active_alphas_sed = {a for a, on in [(1.0002, show_alpha_sed_0002), (1.0003, show_alpha_sed_0003),
                                          (1.0004, show_alpha_sed_0004)] if on}
    # ── α_r active sets (reduction: unmodified 'src' / modified 'sed') ──
    # A custom α_r (if set) is always shown on both protoliths.
    _custom_set = {custom_ar} if custom_ar else set()
    _active_ar_src = {a for a, on in [(1.00002, show_ar_src_0001), (1.0002, show_ar_src_0002),
                                      (1.0003, show_ar_src_0003), (1.0004, show_ar_src_0004)]
                      if on} | _custom_set
    _active_ar_sed = {a for a, on in [(1.00002, show_ar_sed_0001), (1.0002, show_ar_sed_0002),
                                      (1.0003, show_ar_sed_0003), (1.0004, show_ar_sed_0004)]
                      if on} | _custom_set

    def _seg_dbulk(seg):
        """Pick the D_bulk row for the selected salinity and apply the per-segment φ."""
        rdf = seg.get('results_df')
        if rdf is None or rdf.empty:
            return None
        row = rdf[rdf['NaCl_wt_pct'] == selected_salinity]
        if row.empty:
            row = rdf.iloc[[0]]
        def g(c):
            return float(row[c].iloc[0]) if c in row.columns and pd.notna(row[c].iloc[0]) else None
        d = {k: g(k) for k in ('MO_Dbulk', 'CE_Dbulk', 'U_Dbulk', 'TH_Dbulk',
                               'W_Dbulk', 'NB_Dbulk', 'LA_Dbulk')}
        # Apply φ. The bracket was already evaluated at the segment's own T_K,
        # so the bracket and φ share the temperature. φ is the per-segment Excel value.
        # Mo: all minerals corrected. W: only CPX & GRT (fixed W_D_RUT = 1.250 uncorrected).
        if apply_phi600 and d['MO_Dbulk'] is not None:
            phi_seg = seg['phi']
            d['MO_Dbulk'] *= phi_seg
            if d['W_Dbulk'] is not None:
                m_rut_frac = (seg.get('mode_rut') or 0.0) / 100.0
                W_rut_contrib = W_D_RUT_FIXED * m_rut_frac
                d['W_Dbulk'] = (d['W_Dbulk'] - W_rut_contrib) * phi_seg + W_rut_contrib
        return d

    def _seg_results(seg, d):
        """Unaltered-source Stage-1 extraction Rayleigh curves for one segment.

        Returns the FULL dict (all α) — α-filtering for display happens at draw
        time, and the tables show every α. Also used as the source of the
        Stage-2 reduction start point (α=ALPHA_EXT_REF branch at F=F_BRANCH)."""
        C0 = seg['C0']
        seg_d0 = seg['d98Mo'] if seg.get('d98Mo') is not None else delta0
        try:
            res = compute_fluid_fractionation(
                C0_MO=C0.get('MO'), C0_Ce=C0.get('Ce'),
                MO_Dbulk=d['MO_Dbulk'], CE_Dbulk=d['CE_Dbulk'],
                C0_U=C0.get('U'), C0_Th=C0.get('Th'),
                U_Dbulk=d['U_Dbulk'], TH_Dbulk=d['TH_Dbulk'],
                C0_W=C0.get('W'), W_Dbulk=d['W_Dbulk'],
                C0_Nb=C0.get('Nb'), C0_La=C0.get('La'),
                NB_Dbulk=d['NB_Dbulk'], LA_Dbulk=d['LA_Dbulk'],
                delta0=seg_d0, alphas=alphas,
            )
        except (ValueError, ZeroDivisionError):
            return None
        return res

    def _seg_results_sed(seg, d):
        """Sediment-modified Rayleigh curves → (results_sed, delta0_sed, C0_sed)."""
        if not (sed_pct_val and carb_hemi_pct_val is not None and sed_pct_val > 0):
            return None, None, None
        C0 = seg['C0']
        seg_d0 = seg['d98Mo'] if seg.get('d98Mo') is not None else delta0
        f_sed = sed_pct_val / 100.0
        f_hemi = carb_hemi_pct_val / 100.0
        mixed_sed = make_mixed_sediment(f_hemi)
        src_Mo = C0.get('MO') or 0.0
        src_Ce = C0.get('Ce') or 0.0
        src_U = C0.get('U') or 0.0
        src_Th = C0.get('Th') or 0.0
        src_W = C0.get('W') or 0.0
        src_Nb = C0.get('Nb') or 0.0
        src_La = C0.get('La') or 0.0
        C0_MO_sed = (1 - f_sed) * src_Mo + f_sed * mixed_sed.get('Mo', 0.0)
        C0_Ce_sed = (1 - f_sed) * src_Ce + f_sed * mixed_sed.get('Ce', 0.0)
        C0_U_sed = (1 - f_sed) * src_U + f_sed * mixed_sed.get('U', 0.0) if src_U else None
        C0_Th_sed = (1 - f_sed) * src_Th + f_sed * mixed_sed.get('Th', 0.0) if src_Th else None
        C0_W_sed = (1 - f_sed) * src_W + f_sed * mixed_sed.get('W', 0.0) if src_W else None
        C0_Nb_sed = (1 - f_sed) * src_Nb + f_sed * mixed_sed.get('Nb', 0.0) if src_Nb else None
        C0_La_sed = (1 - f_sed) * src_La + f_sed * mixed_sed.get('La', 0.0) if src_La else None
        d98Mo_sed_mix = (
            D98MO_CARBONATES * SUBDUCTED_CARBONATES['Mo'] * (1 - f_hemi)
            + D98MO_HEMIPELAGIC * SUBDUCTED_HEMIPELAGIC['Mo'] * f_hemi
        ) / mixed_sed.get('Mo', 1.0)
        d0_sed = (((1 - f_sed) * src_Mo * seg_d0 + f_sed * mixed_sed.get('Mo', 0.0) * d98Mo_sed_mix)
                  / C0_MO_sed) if C0_MO_sed > 0 else seg_d0
        try:
            res_sed = compute_fluid_fractionation(
                C0_MO=C0_MO_sed, C0_Ce=C0_Ce_sed,
                MO_Dbulk=d['MO_Dbulk'], CE_Dbulk=d['CE_Dbulk'],
                C0_U=C0_U_sed, C0_Th=C0_Th_sed,
                U_Dbulk=d['U_Dbulk'], TH_Dbulk=d['TH_Dbulk'],
                C0_W=C0_W_sed, W_Dbulk=d['W_Dbulk'],
                C0_Nb=C0_Nb_sed, C0_La=C0_La_sed,
                NB_Dbulk=d['NB_Dbulk'], LA_Dbulk=d['LA_Dbulk'],
                delta0=d0_sed, alphas=alphas,
            )
        except (ValueError, ZeroDivisionError):
            return None, None, None
        C0_sed = {'MO': C0_MO_sed, 'Ce': C0_Ce_sed, 'W': C0_W_sed, 'U': C0_U_sed,
                  'Th': C0_Th_sed, 'Nb': C0_Nb_sed, 'La': C0_La_sed}
        return res_sed, d0_sed, C0_sed

    def _seg_reduction(res_stage1, seg, d):
        """Stage-2 redox-filter curves branching off `res_stage1` at F=F_BRANCH.

        ΔFMQ_filter is the filter's ABSOLUTE redox state (a petrological buffer
        floor). Blank/NaN ΔFMQ_filter ⇒ no filter ⇒ Stage 2 OFF for this segment.

        Mechanism: fluid–rock re-equilibration — exit Mo is rock-buffered,
        C_eq = C_rockF / D_bulk(T_seg, ΔFMQ_filter). The filter D is the
        segment's (φ-corrected) slab D rescaled analytically by the exact Bali
        fO₂ exponents: D_Mo×10^{−0.44Δ}; W silicate part ×10^{−0.07Δ} with the
        fixed rutile term preserved. Valid for dilute (undersaturated) fluids;
        a saturation-limited law does not apply to this dilute fluid.
        Returns the FULL dict (all α_r) or None."""
        dfmq_filter = seg.get('dFMQ_filter')
        dfmq_slab = seg.get('dFMQ')
        if (dfmq_filter is None or dfmq_slab is None
                or not np.isfinite(dfmq_filter) or not np.isfinite(dfmq_slab)):
            return None
        dlogfO2_end = float(dfmq_filter) - float(dfmq_slab)
        start = extract_reduction_start(res_stage1, F_branch=F_BRANCH,
                                        alpha_ref=ALPHA_EXT_REF)
        if start is None:
            return None
        try:
            if d is None or d.get('MO_Dbulk') in (None, 0):
                return None
            # rescale the slab-condition D to the filter's absolute fO₂
            D_Mo_f = d['MO_Dbulk'] * 10.0 ** (-EXP_MO * dlogfO2_end)
            D_W_f = None
            if d.get('W_Dbulk') is not None:
                m_rut_frac = (seg.get('mode_rut') or 0.0) / 100.0
                W_rut = W_D_RUT_FIXED * m_rut_frac
                D_W_f = ((d['W_Dbulk'] - W_rut)
                         * 10.0 ** (-EXP_W * dlogfO2_end) + W_rut)
            # Filter-rock Mo: per-segment Excel column 'C0_MO_filter' when
            # present; blank cells fall back to the dashboard input. W is
            # tied by the hemipelagic ratio 1.22/3.96 in both cases.
            _seg_rockF = seg.get('C0_MO_filter')
            if _seg_rockF is not None and np.isfinite(_seg_rockF) and _seg_rockF > 0:
                _rockF_Mo, _rockF_W = float(_seg_rockF), 1.22 * float(_seg_rockF) / 3.96
            else:
                _rockF_Mo, _rockF_W = rockF_Mo, rockF_W
            red = compute_filter_reequilibration(
                start, alpha_r_list=alpha_r_active_list,
                D_Mo_filter=D_Mo_f, D_W_filter=D_W_f,
                C_rockF_Mo=_rockF_Mo, C_rockF_W=_rockF_W,
                dlogfO2_end=dlogfO2_end, iso_mode=iso_mode)
        except (ValueError, ZeroDivisionError):
            return None
        return red or None

    # ── Build per-segment plotting payloads (skip hidden / not-ready segments) ──
    seg_plot = []
    for _i, seg in enumerate(p2_segments):
        _seg_id = seg.get('row') if seg.get('row') is not None else _i
        if _seg_id not in _visible:
            continue
        _phi_src = "from Excel" if seg.get('phi_from_excel') else "computed from T_K"
        if not seg.get('ready'):
            st.info(
                f"**{seg['name']}** — φ = {seg['phi']:.4f} ({_phi_src}), "
                f"T_K = {seg['T_K']:.2f} K. Fill P_bar, ΔFMQ, modal and C0 columns in the "
                f"'Part2_Segments' sheet to plot this segment."
            )
            continue
        d = _seg_dbulk(seg)
        if d is None or d['MO_Dbulk'] is None or d['CE_Dbulk'] is None:
            continue
        res = _seg_results(seg, d)
        if res is None:
            continue
        res_sed, d0_sed, C0_sed = _seg_results_sed(seg, d)
        # Stage-2 redox filter (branches off each Stage-1 extraction at F=F_BRANCH;
        # driven by the segment's ΔFMQ_filter − ΔFMQ_slab; NaN filter ⇒ off)
        res_red = _seg_reduction(res, seg, d)
        res_sed_red = _seg_reduction(res_sed, seg, d) if res_sed is not None else None
        # No-op filter warning: the clip C_eq = min(C_rockF/D, Mo_start) fired at
        # g_eq = 1 → this ΔFMQ_filter removes NO Mo (curve collapses to a point).
        if res_red:
            _any_df = next(iter(res_red.values()))
            _ep = _any_df[_any_df['is_endpoint'] == 1]
            if not _ep.empty and _ep['g'].iloc[0] >= 1.0 - 1e-9:
                st.warning(
                    f"**{seg['name']}** — the redox filter removes NO Mo: at "
                    f"ΔFMQ_filter = {seg.get('dFMQ_filter'):+.2f} (slab "
                    f"ΔFMQ = {seg.get('dFMQ'):+.2f}) the rock-buffered value "
                    f"C_rockF/D_bulk is ≥ the fluid's starting Mo, so the clip "
                    f"(filter never ADDS Mo) makes Stage 2 a no-op. Use a more "
                    f"reducing ΔFMQ_filter or a lower filter-rock Mo."
                )
        seg_plot.append({
            'name': seg['name'],
            'row': seg.get('row'),
            'protolith': seg.get('protolith', seg['name']),
            'color': {'Nicaragua': '#0072B2', 'Costa Rica': '#D55E00'}.get(
                _seg_region(seg['name']), _SEGMENT_COLORS[_i % len(_SEGMENT_COLORS)]),
            'phi': seg['phi'], 'phi_src': _phi_src,
            'P_bar': seg.get('P_bar'), 'T_K': seg['T_K'],
            'dfmq': seg.get('dFMQ'), 'dfmq_filter': seg.get('dFMQ_filter'),
            'mode_rut': seg.get('mode_rut'), 'C0_MO_filter': seg.get('C0_MO_filter'),
            'phi_applied': bool(apply_phi600), 'salinity': selected_salinity,
            'd98Mo': seg.get('d98Mo'), 'd98Mo_2SD': seg.get('d98Mo_2SD'),
            'C0': seg['C0'], 'dbulk': d,
            'results': res, 'results_sed': res_sed, 'delta0_sed': d0_sed,
            'C0_sed': C0_sed,
            'results_red': res_red, 'results_sed_red': res_sed_red,
        })

    if not seg_plot:
        st.info(
            "No ready segments to plot yet. Fill P_bar, ΔFMQ, modal and C0 in the "
            "'Part2_Segments' sheet (φ and T_K are already set)."
        )
        return

    # Active-source readout: protolith, T_K, φ (bracket & φ share T_K), starting δ⁹⁸/⁹⁵Mo.
    def _seg_caption(s):
        _d0, _sd = s.get('d98Mo'), s.get('d98Mo_2SD')
        if _d0 is None:
            _d0txt = ""
        elif _sd is not None:
            _d0txt = f", δ⁹⁸/⁹⁵Mo₀ = {_d0:+.2f} ± {_sd:.2f}‰"
        else:
            _d0txt = f", δ⁹⁸/⁹⁵Mo₀ = {_d0:+.2f}‰"
        return (f"**{s['name']}** (row {s.get('row')}, {s.get('protolith','')}): "
                f"T_K = {s['T_K']:.2f} K, φ = {s['phi']:.4f} ({s['phi_src']}){_d0txt}")
    st.caption("  ·  ".join(_seg_caption(s) for s in seg_plot) + "  ·  φ ref = 1000 °C")

    # The ΔFMQ/rutile annotation uses the first segment's values (segments share the plot box).
    dfmq_val = seg_plot[0]['dfmq']
    mode_rut_pct = seg_plot[0]['mode_rut']

    # ========== 3×2 GRID OF PLOTS ==========
    # Toolbar image-export settings
    _p2_dl_config = {
        'toImageButtonOptions': {
            'format': 'png',
            'filename': 'bali_part2_plot',
            'height': 650,
            'width': 560,
            'scale': 4
        },
        'displayModeBar': True,
        'displaylogo': False
    }

    def _draw_segment_curves(fig, x_col):
        """Overlay the active source's TWO-STAGE, TWO-PROTOLITH Rayleigh curve system.

        Colour = stage (Extraction → red, Reduction → blue).
        Line    = protolith (Unmodified → solid, Modified → dashed).
        α       = one trace each, α shown in legend + hover.
        Extraction curves carry F% markers (10/20/30); reduction curves carry the
        removal markers and branch off the extraction curve at F_BRANCH.

        Visibility respects the protolith toggles (show_unmodified / show_modified),
        the stage toggles (show_extraction / show_reduction) and the per-protolith
        active-α sets (extraction: _active_alphas_src/_sed; reduction: _active_ar_src/_sed).
        """
        _hemi = int(round(carb_hemi_pct_val)); _carb = 100 - _hemi; _sed = int(round(sed_pct_val))

        def _source_x(c0):
            """Protolith source x-coordinate in this plot's x_col space."""
            c0 = c0 or {}
            def rt(a, b):
                va, vb = c0.get(a), c0.get(b)
                return (va / vb) if (va and vb) else None
            xmap = {'cf_cum_Mo': c0.get('MO'),
                    'MoCe_f_cum': rt('MO', 'Ce'), 'CeMo_f_cum': rt('Ce', 'MO'),
                    'UTh_f_cum': rt('U', 'Th'), 'WTh_f_cum': rt('W', 'Th'),
                    'MoW_f_cum': rt('MO', 'W'), 'WCe_f_cum': rt('W', 'Ce'),
                    'NbLa_f_cum': rt('Nb', 'La')}
            return xmap.get(x_col)

        def _add_source_star(x, y, label, with_x):
            """Yellow source star (thin black border). The MODIFIED protolith
            carries a small black × overlaid at the centre to distinguish it."""
            if x is None or not np.isfinite(x) or x <= 0 or y is None:
                return
            fig.add_trace(go.Scatter(
                x=[x], y=[y], mode='markers', name=label,
                marker=dict(color='#FFD700', size=13, symbol='star',
                            line=dict(color='black', width=0.8)),
                hovertemplate=(f'{label}<br>x = %{{x:.4f}}<br>'
                               'δ⁹⁸/⁹⁵Mo = %{y:.3f}‰<extra></extra>'),
            ))
            if with_x:
                fig.add_trace(go.Scatter(
                    x=[x], y=[y], mode='markers', showlegend=False,
                    name=f'{label} ×',
                    marker=dict(color='black', size=4, symbol='x-thin',
                                line=dict(color='black', width=1.2)),
                    hoverinfo='skip',
                ))

        for s in seg_plot:
            _ab = ("Nic" if "nicaragua" in s['name'].lower()
                   else ("CR" if "costa" in s['name'].lower() else s['name']))
            _sed_tag = f"{_sed}%sed {_carb}C:{_hemi}H"
            has_sed = s['results_sed'] is not None

            # ── Protolith source starts (yellow stars) ──
            # Unmodified: plain star. Modified (sediment-mixed): star + black ×.
            if show_unmodified:
                _d0u = s.get('d98Mo') if s.get('d98Mo') is not None else -0.21
                _add_source_star(_source_x(s.get('C0')), _d0u,
                                 f"{_ab} protolith source", with_x=False)
            if show_modified and has_sed and s.get('C0_sed'):
                _add_source_star(_source_x(s.get('C0_sed')), s.get('delta0_sed'),
                                 f"{_ab} modified source ({_sed_tag})", with_x=True)

            # ── Unmodified protolith (solid) ──
            if show_unmodified:
                if show_extraction and _active_alphas_src:
                    _add_stage_curves(
                        fig, {a: df for a, df in s['results'].items() if a in _active_alphas_src},
                        x_col, _STAGE_COLOR_EXTRACT, 'solid',
                        f"{_ab} UndM · Extr", 'F', F_HIGHLIGHT)
                if show_reduction and _active_ar_src and s['results_red']:
                    _add_stage_curves(
                        fig, {a: df for a, df in s['results_red'].items() if a in _active_ar_src},
                        x_col, _STAGE_COLOR_REDUCE, 'solid',
                        f"{_ab} UndM · Redn[{'batch' if iso_mode == 'batch' else 'Rayl'}]",
                        'g', G_HIGHLIGHT)

            # ── Sediment-modified protolith (dashed) ──
            if show_modified and has_sed:
                if show_extraction and _active_alphas_sed:
                    _add_stage_curves(
                        fig, {a: df for a, df in s['results_sed'].items() if a in _active_alphas_sed},
                        x_col, _STAGE_COLOR_EXTRACT, 'dash',
                        f"{_ab} Mod ({_sed_tag}) · Extr", 'F', F_HIGHLIGHT)
                if show_reduction and _active_ar_sed and s['results_sed_red']:
                    _add_stage_curves(
                        fig, {a: df for a, df in s['results_sed_red'].items() if a in _active_ar_sed},
                        x_col, _STAGE_COLOR_REDUCE, 'dash',
                        f"{_ab} Mod ({_sed_tag}) · Redn[{'batch' if iso_mode == 'batch' else 'Rayl'}]",
                        'g', G_HIGHLIGHT)
        return fig
    
    # ========== PLOTS ==========
    # Default: only Plot D (δ⁹⁸/⁹⁵Mo vs Mo/W) and Plot H (δ⁹⁸/⁹⁵Mo vs Mo) are
    # computed and rendered. The other plots are OFF by default and are NOT
    # computed while hidden — this makes parameter changes much faster.
    _s0 = seg_plot[0]
    _fil = _s0.get('dfmq_filter')
    _info_lines = [
        f"{_s0['name']} slab segment, N-MORB source",
        f"ΔFMQ<sub>slab</sub> = {_s0['dfmq']:+g} ({_s0['T_K']:.2f} K)",
        (f"ΔFMQ<sub>filter</sub> = {_fil:+g} (f<sub>branch</sub> = {F_BRANCH:g})"
         if _fil is not None and np.isfinite(_fil) else "no redox filter"),
        f"{salinity:g} wt.% NaCl",
        f"{mode_rut_pct:g}% Rut.",
    ]
    with _plot_col:
        st.markdown("**δ⁹⁸/⁹⁵Mo of the cumulative slab fluid against its Mo content**")
        fig_mo = create_fractionation_plot_vs_mo(
            title="<br>".join(_info_lines),
            log_x=True,
            arc_df=arc_df,
            style_df=style_df,
            arc_settings=arc_settings,
            dfmq=dfmq_val,
            mode_rut_pct=mode_rut_pct
        )
        fig_mo = _draw_segment_curves(fig_mo, 'cf_cum_Mo')
        st.plotly_chart(fig_mo, use_container_width=False, config=_p2_dl_config, key="p2_plot_mo")

    if st.button("Prepare figure downloads (PNG 300 dpi + SVG vector)",
                 key="p2_prepare_figure_downloads"):
        st.caption("Export plots as PNG (1500×1500 px ≈ 300 dpi at 5 inch) or SVG (scalable vector).")
        
        # Main plots (always available)
        st.markdown("**δ⁹⁸/⁹⁵Mo vs Mo**")
        dl_col8a, dl_col8b = st.columns(2)

        with dl_col8a:
            try:
                png_mo = fig_mo.to_image(format="png", width=1500, height=1500, scale=1)
                st.download_button(
                    label="PNG (Mo)",
                    data=png_mo,
                    file_name="two_stage_model_d98Mo_vs_Mo.png",
                    mime="image/png",
                    key="dl_mo_png",
                    use_container_width=True
                )
            except Exception as e:
                st.error(f"PNG export error: {str(e)[:100]}")

        with dl_col8b:
            try:
                svg_mo = fig_mo.to_image(format="svg", width=560, height=650)
                st.download_button(
                    label="SVG (Mo)",
                    data=svg_mo,
                    file_name="two_stage_model_d98Mo_vs_Mo.svg",
                    mime="image/svg+xml",
                    key="dl_mo_svg",
                    use_container_width=True
                )
            except Exception as e:
                st.error(f"SVG export error: {str(e)[:100]}")

    # ── Parameters and tables of the active segment ──
    _rep = seg_plot[0]   # active slab segment
    delta0 = _rep['d98Mo'] if _rep.get('d98Mo') is not None else delta0

    _line_table = build_plot_line_table(fig_mo)
    if not _line_table.empty:
        with st.expander("Copy exact plotted X/Y pairs of the displayed lines (Excel / DataGraph)"):
            st.caption("X = Mo (µg/g), Y = δ⁹⁸/⁹⁵Mo (‰); one X/Y pair per displayed curve, points in plotting order. "
                       "Use the Copy button of the code block to keep headers and precision.")
            _line_tsv = _line_table.to_csv(sep='\t', index=False, na_rep='', float_format='%.17g',
                                           lineterminator='\n')
            st.code(_line_tsv, language=None)
            st.download_button("Download as TSV", data=_line_tsv, mime="text/tab-separated-values",
                               file_name=f"two_stage_model_row{_rep.get('row')}_lines.tsv",
                               key="p2_datagraph_lines_download")

    # ========== DETAILED CALCULATION TABLES ==========
    # Four tables for the active source: {Unmodified, Modified} × {Extraction,
    # Reduction}. Each table has ONE column per α (all α included, even hidden
    # ones). Recomputed on every slider / checkbox change.
    st.markdown("#### Detailed Calculation Tables")
    st.caption("Active source: "
               f"**{_rep['name']}** (row {_rep.get('row')}, {_rep.get('protolith','')}). "
               "Unmodified & sediment-modified protoliths × extraction / reduction stages, "
               "one δ⁹⁸/⁹⁵Mo column per α.")

    _t_hemi = int(round(carb_hemi_pct_val)); _t_carb = 100 - _t_hemi; _t_sed = int(round(sed_pct_val))
    _sed_tag = f"{_t_sed}%sed, {_t_carb}C:{_t_hemi}H"

    def _stage_table(res_dict, progress_col):
        """Compact per-stage table: rows = F (extraction) or g (reduction);
        α-independent shared columns once; one δ column per α (all α)."""
        if not res_dict:
            return None
        _as = sorted(res_dict.keys())
        base = res_dict[_as[0]]
        if progress_col not in base.columns:
            return None
        out = pd.DataFrame()
        out[progress_col] = base[progress_col].values
        if progress_col == 'g':
            # Redox-filter stage: g is EMERGENT from the fO₂ drop. Show the
            # equivalent Δlog fO₂, the removal %, the W survival, and the
            # precipitate mass-balance ledger (fluid + precipitate = start).
            if 'fex' in base.columns:
                out['f_ex'] = base['fex'].values
            if 'dlogfO2' in base.columns:
                out['Δlog fO₂'] = base['dlogfO2'].values
            if 'removal_pct' in base.columns:
                out['reduced_%'] = base['removal_pct'].values
            else:
                out['reduced_%'] = (1.0 - base['g'].astype(float)) * 100.0
            _shared = [
                'Mo_start', 'Mo_eq_unclipped', 'Mo_eq', 'D_Mo_filter',
                'C_rockF_Mo', 'delta_start', 'W_start', 'W_eq_unclipped',
                'W_eq', 'D_W_filter', 'C_rockF_W',
                'cf_cum_Mo', 'Mo_precip_cum', 'cf_cum_Ce', 'cf_cum_U',
                'cf_cum_Th', 'cf_cum_W', 'cf_cum_Nb', 'cf_cum_La',
                'g_W', 'MoW_f_cum', 'MoCe_f_cum', 'UTh_f_cum',
                'WTh_f_cum', 'NbLa_f_cum', 'WCe_f_cum', 'CeMo_f_cum',
            ]
        else:
            _shared = [
                'cs_Mo', 'cf_cum_Mo', 'fs', 'cf_cum_Ce', 'cf_cum_U',
                'cf_cum_Th', 'cf_cum_W', 'cf_cum_Nb', 'cf_cum_La',
                'MoCe_f_cum', 'CeMo_f_cum', 'UTh_f_cum', 'WTh_f_cum',
                'MoW_f_cum', 'NbLa_f_cum', 'WCe_f_cum',
            ]
        for c in _shared:
            if c in base.columns:
                out[c] = base[c].values
        for a in _as:
            out[f'δ(α={a})'] = res_dict[a]['delta_f_cum'].values
        if progress_col == 'g':
            for a in _as:
                if 'delta_precip_cum' in res_dict[a].columns:
                    out[f'δ_precip(α={a})'] = res_dict[a]['delta_precip_cum'].values
        return out

    _tbl_groups = [
        ("Unmodified · Extraction (Stage 1)", _rep.get('results'), 'F'),
        ("Unmodified · Reduction (Stage 2)", _rep.get('results_red'), 'g'),
    ]
    if _rep.get('results_sed') is not None:
        _tbl_groups.append((f"Modified ({_sed_tag}) · Extraction (Stage 1)",
                            _rep.get('results_sed'), 'F'))
    if _rep.get('results_sed_red') is not None:
        _tbl_groups.append((f"Modified ({_sed_tag}) · Reduction (Stage 2)",
                            _rep.get('results_sed_red'), 'g'))

    for _gi, (_gname, _rdict, _pcol) in enumerate(_tbl_groups):
        _tbl = _stage_table(_rdict, _pcol)
        with st.expander(f" {_gname}", expanded=(_gi == 0)):
            if _tbl is None or _tbl.empty:
                st.info("No data for this protolith × stage at the current settings.")
                continue
            _stage_kind = (("Reduction (Stage 2, redox filter — mass: C_eq = "
                            "C_rockF/D(T, ΔFMQ_filter); isotopes ["
                            + ("BATCH closed-system: δ = δ_start + (1−g)·1000·lnα_r"
                               if iso_mode == "batch" else
                               "RAYLEIGH open-system: δ = δ_start − 1000·lnα_r·ln g")
                            + "] — removed Mo is LIGHT, fluid drifts HEAVIER)")) if _pcol == 'g' \
                else "Extraction (Stage 1, fluid gains Mo)"
            _alpha_lbl = "α_r" if _pcol == 'g' else "α"
            st.caption(f"{_stage_kind} — {_alpha_lbl} ∈ "
                       f"{{{', '.join(str(a) for a in sorted(_rdict.keys()))}}}. "
                       f"δ⁹⁸/⁹⁵Mo in ‰; one column per {_alpha_lbl}. "
                       + ("g = Mo fraction remaining in the fluid, EMERGENT from "
                          "re-equilibration with the filter rock (C_eq = "
                          "C_rockF/D(T, ΔFMQ_filter), clipped ≤ start); branches off "
                          f"extraction at F = {F_BRANCH:.0%} "
                          f"(extraction α = {ALPHA_EXT_REF}). "
                          "δ_precip = stay-behind reduced "
                          "precipitate (mass balance: fluid + precipitate = start)."
                          if _pcol == 'g'
                          else "F = cumulative fluid mass fraction extracted."))
            _fmt = {}
            for c in _tbl.columns:
                if c in ('F', 'g'):
                    _fmt[c] = '{:.2f}'
                elif c == 'f_ex':
                    _fmt[c] = '{:.2f}'
                elif c == 'reduced_%':
                    _fmt[c] = '{:.0f}'
                elif c == 'Δlog fO₂':
                    _fmt[c] = '{:.3f}'
                elif c == 'alpha_r':
                    _fmt[c] = '{:.7f}'
                elif c == 'isotope_mode':
                    continue
                elif c.startswith('δ'):
                    _fmt[c] = '{:.4f}'
                elif c in ('fs', 'cf_cum_Mo', 'cs_Mo', 'g_W', 'Mo_precip_cum',
                           'Mo_start', 'Mo_eq_unclipped', 'Mo_eq', 'C_rockF_Mo',
                           'W_start', 'W_eq_unclipped', 'W_eq', 'C_rockF_W',
                           'cf_cum_Ce', 'cf_cum_U', 'cf_cum_Th', 'cf_cum_W',
                           'cf_cum_Nb', 'cf_cum_La'):
                    _fmt[c] = '{:.6f}'
                else:
                    _fmt[c] = '{:.5f}'
            st.dataframe(_tbl.style.format(_fmt, na_rep='NaN'),
                         use_container_width=True, height=360)
            if _pcol == 'g':
                _src_C0 = _rep.get('C0') or {}
                _used_C0 = (_rep.get('C0_sed') if 'Modified' in _gname else _src_C0) or {}
                _src_delta = (_rep.get('delta0_sed') if 'Modified' in _gname
                              else _rep.get('d98Mo'))
                _dbulk = _rep.get('dbulk') or {}
                _context = {
                    'segment': _rep.get('name'),
                    'source_row': _rep.get('row'),
                    'protolith': _rep.get('protolith'),
                    'P_bar': _rep.get('P_bar'),
                    'T_K': _rep.get('T_K'),
                    'phi': _rep.get('phi'),
                    'phi_applied': _rep.get('phi_applied'),
                    'mode_rut_pct': _rep.get('mode_rut'),
                    'NaCl_wt_pct': _rep.get('salinity'),
                    'DeltaFMQ_slab': _rep.get('dfmq'),
                    'DeltaFMQ_filter': _rep.get('dfmq_filter'),
                    'C0_MO_used': _used_C0.get('MO'),
                    'C0_Ce_used': _used_C0.get('Ce'),
                    'C0_U_used': _used_C0.get('U'),
                    'C0_Th_used': _used_C0.get('Th'),
                    'C0_W_used': _used_C0.get('W'),
                    'C0_Nb_used': _used_C0.get('Nb'),
                    'C0_La_used': _used_C0.get('La'),
                    'delta98Mo_source_used': _src_delta,
                    'D_Mo_slab_phi_corrected': _dbulk.get('MO_Dbulk'),
                    'D_Ce_slab_phi_corrected': _dbulk.get('CE_Dbulk'),
                    'D_U_slab_phi_corrected': _dbulk.get('U_Dbulk'),
                    'D_Th_slab_phi_corrected': _dbulk.get('TH_Dbulk'),
                    'D_W_slab_phi_corrected': _dbulk.get('W_Dbulk'),
                    'D_Nb_slab_phi_corrected': _dbulk.get('NB_Dbulk'),
                    'D_La_slab_phi_corrected': _dbulk.get('LA_Dbulk'),
                    'F_branch': F_BRANCH,
                    'alpha_extraction_reference': ALPHA_EXT_REF,
                    'isotope_mode': iso_mode,
                    'C0_MO_filter_segment': _rep.get('C0_MO_filter'),
                    'sediment_pct': _t_sed if 'Modified' in _gname else 0,
                    'carb_hemi_pct': _t_hemi if 'Modified' in _gname else np.nan,
                }
                _export_tbl = _tbl.assign(**_context)
                _tsv = _export_tbl.to_csv(
                    sep='\t', index=False, na_rep='NaN', float_format='%.12g',
                    lineterminator='\n')
                _base_filename = f"part2_stage2_row{_rep.get('row')}_table{_gi + 1}"
                st.download_button(
                    "Download tab-separated data (Excel / DataGraph)",
                    data=_tsv,
                    file_name=f"{_base_filename}.tsv",
                    mime="text/tab-separated-values",
                    key=f"p2_stage2_tsv_{_rep.get('row')}_{_gi}",
                )
                with st.expander("Copy/paste tab-separated values with model parameters"):
                    st.caption(
                        "Copy the tab-delimited text from the code block, or use the "
                        "download button above. It contains full-precision model outputs "
                        "and the fixed inputs repeated on every row."
                    )
                    st.code(_tsv, language=None)
