"""
Fluid–mineral partitioning engine of Bali et al. (2011, 2012), used by the
two-stage fractionation model to obtain bulk partition coefficients for each
slab segment at its own temperature and oxygen fugacity.

Contents: model constants, normalisation of the input-workbook headers, the
FMQ buffer, and the vectorised partitioning calculation (one row per input row
and salinity).
"""

import math

import numpy as np
import pandas as pd


class BaliConstants:
    """Constants of the partitioning model (concentrations in µg/g)."""

    def __init__(self):
        self.molar_masses = {"MO": 95.95, "W": 183.84, "NACL": 58.44}
        # Depleted mantle and N-MORB, Table B1 of the article
        self.dm_concentrations = {
            "MO": 0.025, "Ce": 0.772, "W": 0.002, "U": 0.005,
            "Th": 0.014, "Nb": 0.21, "La": 0.234
        }
        self.morb_concentrations = {
            "MO": 0.408, "Ce": 12.0, "W": 0.038, "U": 0.071,
            "Th": 0.187, "Nb": 3.507, "La": 3.895
        }
        # log10 of the molality → µg/g scale factors
        self.log_scales = {
            "MO": math.log10(self.molar_masses["MO"] * 1000.0),
            "W": math.log10(self.molar_masses["W"] * 1000.0)
        }


CONSTANTS = BaliConstants()

# Module-level aliases used by the vectorised calculation
DM = CONSTANTS.dm_concentrations
M_MO = CONSTANTS.molar_masses["MO"]
M_W = CONSTANTS.molar_masses["W"]
M_NACL = CONSTANTS.molar_masses["NACL"]
CE_D_CPX, CE_D_GRT, CE_D_RUT = 2.0, 0.4, 2.0
TH_D_CPX, TH_D_GRT, TH_D_RUT = 1.19, 0.610, 0.10
NB_D_CPX, NB_D_GRT, NB_D_RUT = 0.172, 0.204, 200.0
LA_D_CPX, LA_D_GRT, LA_D_RUT = 1.429, 0.204, 1.250
W_D_RUT = 1.250
SCALE_LOG10_M_MO = CONSTANTS.log_scales["MO"]
SCALE_LOG10_M_W = CONSTANTS.log_scales["W"]

# Header normalisation of the input workbook
def _norm(s: str) -> str:
    if not isinstance(s, str): return ""
    x = s.strip().lower()
    # Normalize common delta characters to ascii 'd' so headers like 'ΔFMQ' map to 'dfmq'
    x = x.replace('Δ', 'd').replace('δ', 'd')
    for ch in ("°", " ", "\t", "(", ")", "_", "-"):
        x = x.replace(ch, "")
    return x

HEADER_MAP = {
    # identifiers
    "input": "rock_id", "id": "rock_id", "sample": "rock_id",
    # pressure (BAR)
    "pbar": "P_bar", "p(bar)": "P_bar", "pressurebar": "P_bar", "p": "P_bar",
    # temperature
    "temperaturek": "T_K", "temperature(k)": "T_K", "tk": "T_K",
    "temperaturec": "T_C", "temperature(c)": "T_C", "tc": "T_C",
    # ΔFMQ
    "Δfmq(logunits)": "dFMQ", "dfmq(logunits)": "dFMQ", "Δfmq": "dFMQ", "dfmq": "dFMQ",
    # modal proportion fallbacks (normalized keys)
    "modegrt": "mode_grt", "grt": "mode_grt", "garnet": "mode_grt",
    "modecpx": "mode_cpx", "cpx": "mode_cpx", "clinopyroxene": "mode_cpx",
    "moderut": "mode_rut", "rut": "mode_rut", "rutile": "mode_rut",
    # simple trace element column mappings (if present as exact names)
    "mo": "C0_MO", "molibdenum": "C0_MO",
    "ce": "C0_Ce", "cerium": "C0_Ce",
    "w": "C0_W", "tungsten": "C0_W",
    "u": "C0_U", "uranium": "C0_U",
    "th": "C0_Th", "thorium": "C0_Th",
    "nb": "C0_Nb", "niobium": "C0_Nb",
    "la": "C0_La", "lanthanum": "C0_La",
}


def load_excel(path: str, sheet=0) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name=sheet, engine='openpyxl')

def normalize_headers(df: pd.DataFrame) -> pd.DataFrame:
    # build normalized key map for original columns
    norm_map = {col: _norm(col) for col in df.columns}

    ren = {}
    for col, key in norm_map.items():
        if key in HEADER_MAP:
            ren[col] = HEADER_MAP[key]
    out = df.rename(columns=ren).copy()

    # Fallback: if dFMQ wasn't matched, look for any column whose normalized name contains 'dfmq' or 'fmq'
    if 'dFMQ' not in out.columns:
        for col, key in norm_map.items():
            if 'dfmq' in key or ('fmq' in key and 'd' in key):
                out['dFMQ'] = df[col]
                break
        else:
            # also accept columns that include 'fmq' anywhere
            for col, key in norm_map.items():
                if 'fmq' in key:
                    out['dFMQ'] = df[col]
                    break

    # Ensure rock_id
    if "rock_id" not in out.columns:
        if "input" in out.columns:
            out["rock_id"] = out["input"].astype(str)
        else:
            out.insert(0, "rock_id", [f"rock_{i+1}" for i in range(len(out))])

    # Celsius -> Kelvin if needed
    if "T_K" not in out.columns and "T_C" in out.columns:
        out["T_K"] = pd.to_numeric(out["T_C"], errors="coerce") + 273.15

    # Try to populate modal and trace columns from flexible header variants using norm_map
    # This keeps backward compatibility with the row-wise find_val approach but creates proper
    # dataframe columns for vectorized processing.
    def _assign_first_variant(variants, colname):
        if colname in out.columns:
            return
        for col, key in norm_map.items():
            if key in variants or any(v in key for v in variants):
                out[colname] = df[col]
                return

    _assign_first_variant({'modegrt', 'grt', 'garnet'}, 'mode_grt')
    _assign_first_variant({'modecpx', 'cpx', 'clinopyroxene'}, 'mode_cpx')
    _assign_first_variant({'moderut', 'rut', 'rutile'}, 'mode_rut')

    _assign_first_variant({'mo', 'molibdenum', 'c0moppm', 'c0_mo_ppm', 'c0_mo'}, 'C0_MO')
    _assign_first_variant({'ce', 'cerium', 'c0ceppm', 'c0_ce_ppm', 'c0_ce'}, 'C0_Ce')
    _assign_first_variant({'w', 'tungsten', 'c0wppm', 'c0_w_ppm', 'c0_w'}, 'C0_W')
    _assign_first_variant({'u', 'uranium', 'c0uppm', 'c0_u_ppm', 'c0_u'}, 'C0_U')
    _assign_first_variant({'th', 'thorium', 'c0thppm', 'c0_th_ppm', 'c0_th'}, 'C0_Th')
    _assign_first_variant({'nb', 'niobium', 'c0nbppm', 'c0_nb_ppm', 'c0_nb'}, 'C0_Nb')
    _assign_first_variant({'la', 'lanthanum', 'c0lappm', 'c0_la_ppm', 'c0_la'}, 'C0_La')

    # Coerce numerics (handle unicode minus and stray text) - robust for negative dFMQ
    for c in ("P_bar", "T_K", "dFMQ", 'mode_grt', 'mode_cpx', 'mode_rut', 'C0_MO', 'C0_Ce', 'C0_W', 'C0_U', 'C0_Th', 'C0_Nb', 'C0_La'):
        if c in out.columns:
            # convert to string safely
            s = out[c].astype(str)
            # replace unicode minus (U+2212) with ASCII hyphen-minus and strip whitespace
            s = s.str.replace('\u2212', '-', regex=False).str.strip()
            # allow decimal commas -> convert to dot
            s = s.str.replace(',', '.', regex=False)
            # common textual NA values -> empty string
            s = s.replace({'none': '', 'na': '', 'n/a': '', 'nan': '', '': ''})
            # finally coerce to numeric (invalid -> NaN)
            out[c] = pd.to_numeric(s, errors='coerce')

    # Required columns
    for req in ("rock_id", "P_bar", "T_K"):
        if req not in out.columns:
            raise ValueError(f"Missing required column: {req}")

    return out


# FMQ buffer
def compute_base_fmq_log10(T_K: float, P_bar: float) -> float:
    """
    FMQ base (log10) with P input in BAR (do NOT convert to kbar):
      base = (-25096.3 / T_K) + 8.735 + 0.11 * ((P_bar - 1.0) / T_K)

    This implementation validates T_K and returns NaN for invalid/non-positive T_K values
    to avoid unexpected exceptions or infinite results.
    """
    try:
        t = float(T_K)
        if not math.isfinite(t) or t <= 0.0:
            return float('nan')
    except Exception:
        return float('nan')
    try:
        p_bar = float(P_bar)
    except Exception:
        p_bar = float('nan')
    return (-25096.3 / t) + 8.735 + 0.11 * ((p_bar - 1.0) / t)


# Vectorized helper for FMQ base values (returns numpy array-like)
def compute_base_fmq_log10_vectorized(T_K_arr, P_bar_arr):
    """
    Vectorized FMQ base calculation. Returns an array of float values with NaN
    where inputs are invalid (non-finite or non-positive T_K).
    """
    TK = np.array(T_K_arr, dtype=float)
    PB = np.array(P_bar_arr, dtype=float)
    out = np.full_like(TK, np.nan, dtype=float)
    mask = np.isfinite(TK) & (TK > 0.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        out[mask] = (-25096.3 / TK[mask]) + 8.735 + 0.11 * ((PB[mask] - 1.0) / TK[mask])
    return out



def process_all_vectorized(df: pd.DataFrame, sal_wts):
    """
    Vectorized (but dtype-safe) processing for all rows and a list of salinity wt.% values.
    Returns a pandas DataFrame with one row per input row per salinity.
    This version builds a per-salinity concatenation to preserve pandas dtypes
    and then coerces numeric columns explicitly to avoid object/string issues.
    """
    if df is None or df.empty:
        return pd.DataFrame()

    # Work on a copy
    base_df = df.copy()

    # Ensure expected columns exist
    expected = ['mode_grt', 'mode_cpx', 'mode_rut', 'C0_MO', 'C0_Ce', 'C0_W', 'C0_U', 'C0_Th', 'C0_Nb', 'C0_La', 'dFMQ']
    for c in expected:
        if c not in base_df.columns:
            base_df[c] = np.nan

    # Compute FMQ base and provisional absolute fO2 (dFMQ may be NaN)
    # Use the shared, validated vectorized helper to ensure consistent behavior
    base_df['log10_fO2_FMQ'] = compute_base_fmq_log10_vectorized(base_df['T_K'], base_df['P_bar'])
    base_df['log10_fO2_abs'] = base_df['log10_fO2_FMQ'] + base_df['dFMQ'].fillna(0.0)
    base_df.loc[base_df['dFMQ'].isna(), 'log10_fO2_abs'] = np.nan

    # Build per-salinity frames and concatenate to preserve dtypes
    frames = []
    for wt in sal_wts:
        tmp = base_df.copy()
        tmp['NaCl_wt_pct'] = float(wt)
        frames.append(tmp)
    rep_df = pd.concat(frames, ignore_index=True)

    # Coerce numeric columns to float explicitly to avoid object dtypes
    num_cols = ['P_bar', 'T_K', 'dFMQ', 'mode_grt', 'mode_cpx', 'mode_rut',
                'C0_MO', 'C0_Ce', 'C0_W', 'C0_U', 'C0_Th', 'C0_Nb', 'C0_La',
                'log10_fO2_FMQ', 'log10_fO2_abs']
    for c in num_cols:
        if c in rep_df.columns:
            rep_df[c] = pd.to_numeric(rep_df[c], errors='coerce')

    # molality and log10_NaCl
    wt_series = rep_df['NaCl_wt_pct'].astype(float)
    rep_df['NaCl_m'] = 1000.0 * wt_series / (float(M_NACL) * (100.0 - wt_series))
    rep_df['NaCl_m'] = rep_df['NaCl_m'].replace([np.inf, -np.inf], np.nan)
    rep_df['log10_NaCl'] = np.where((rep_df['NaCl_m'] > 0) & np.isfinite(rep_df['NaCl_m']), np.log10(rep_df['NaCl_m']), np.nan)

    # LogMO
    rep_df['LogMO'] = (0.44 * rep_df['log10_fO2_abs']) + (0.42 * rep_df['log10_NaCl']) - (1.8 * (1000.0 / rep_df['T_K'])) + 4.8

    # E_MO with overflow guards
    scale_log10 = float(SCALE_LOG10_M_MO) if np.isfinite(SCALE_LOG10_M_MO) else np.nan
    combined_log = rep_df['LogMO'] + scale_log10
    ok = np.isfinite(combined_log) & (combined_log <= 308) & (combined_log >= -308)
    rep_df['E_MO'] = np.nan
    rep_df.loc[ok, 'E_MO'] = (10.0 ** rep_df.loc[ok, 'LogMO']) * M_MO * 1000.0
    rep_df['MO_D_CPX'] = np.where(rep_df['E_MO'] > 0, 40.0 / rep_df['E_MO'], np.nan)
    rep_df['MO_D_GRT'] = np.where(rep_df['E_MO'] > 0, 12.0 / rep_df['E_MO'], np.nan)
    rep_df['MO_D_RUT'] = np.where(rep_df['E_MO'] > 0, 87670.0 / rep_df['E_MO'], np.nan)

    # LogW and W partitioning
    rep_df['LogW'] = (0.07 * rep_df['log10_fO2_abs']) - (4.7236 * (1000.0 / rep_df['T_K'])) + 4.4271
    scale_log10_w = float(SCALE_LOG10_M_W) if np.isfinite(SCALE_LOG10_M_W) else np.nan
    combined_log_w = rep_df['LogW'] + scale_log10_w
    okw = np.isfinite(combined_log_w) & (combined_log_w <= 308) & (combined_log_w >= -308)
    rep_df['E_W'] = np.nan
    rep_df.loc[okw, 'E_W'] = (10.0 ** rep_df.loc[okw, 'LogW']) * M_W * 1000.0
    rep_df['W_D_CPX'] = np.where(rep_df['E_W'] > 0, 60.0 / rep_df['E_W'], np.nan)
    rep_df['W_D_GRT'] = np.where(rep_df['E_W'] > 0, 12.0 / rep_df['E_W'], np.nan)

    # LogU & U partitioning
    rep_df['LogU_cg'] = 2.681 + (0.1433 * rep_df['log10_fO2_abs']) + (0.594 * rep_df['NaCl_m'])
    rep_df['LogU_rut'] = 1.7954 + (0.1433 * rep_df['log10_fO2_abs']) + (0.594 * rep_df['NaCl_m'])
    rep_df['E_U_cg'] = np.where(np.isfinite(rep_df['LogU_cg']) & (rep_df['LogU_cg'] <= 308) & (rep_df['LogU_cg'] >= -308), 10.0 ** rep_df['LogU_cg'], np.nan)
    rep_df['E_U_rut'] = np.where(np.isfinite(rep_df['LogU_rut']) & (rep_df['LogU_rut'] <= 308) & (rep_df['LogU_rut'] >= -308), 10.0 ** rep_df['LogU_rut'], np.nan)
    rep_df['U_D_CPX'] = np.where(rep_df['E_U_cg'] > 0, 11.0 / rep_df['E_U_cg'], np.nan)
    rep_df['U_D_GRT'] = np.where(rep_df['E_U_cg'] > 0, 40.0 / rep_df['E_U_cg'], np.nan)
    rep_df['U_D_RUT'] = np.where(rep_df['E_U_rut'] > 0, 94.0 / rep_df['E_U_rut'], np.nan)

    # Compute bulk coefficient (modal-weighted)
    for m in ('mode_cpx', 'mode_grt', 'mode_rut'):
        if m in rep_df.columns:
            rep_df[m] = rep_df[m].astype(float).fillna(0.0)
    m_cpx = rep_df.get('mode_cpx', pd.Series(0.0, index=rep_df.index))
    m_grt = rep_df.get('mode_grt', pd.Series(0.0, index=rep_df.index))
    m_rut = rep_df.get('mode_rut', pd.Series(0.0, index=rep_df.index))

    rep_df['MO_Dbulk'] = (rep_df['MO_D_CPX'].fillna(0.0) * m_cpx) + (rep_df['MO_D_GRT'].fillna(0.0) * m_grt) + (rep_df['MO_D_RUT'].fillna(0.0) * m_rut)
    rep_df['CE_Dbulk'] = (CE_D_CPX * m_cpx) + (CE_D_GRT * m_grt) + (CE_D_RUT * m_rut)
    rep_df['W_Dbulk'] = (rep_df.get('W_D_CPX', 0.0).fillna(0.0) * m_cpx) + (rep_df.get('W_D_GRT', 0.0).fillna(0.0) * m_grt) + (W_D_RUT * m_rut)
    rep_df['U_Dbulk'] = (rep_df.get('U_D_CPX', 0.0).fillna(0.0) * m_cpx) + (rep_df.get('U_D_GRT', 0.0).fillna(0.0) * m_grt) + (rep_df.get('U_D_RUT', 0.0).fillna(0.0) * m_rut)
    rep_df['TH_Dbulk'] = (TH_D_CPX * m_cpx) + (TH_D_GRT * m_grt) + (TH_D_RUT * m_rut)
    rep_df['NB_Dbulk'] = (NB_D_CPX * m_cpx) + (NB_D_GRT * m_grt) + (NB_D_RUT * m_rut)
    rep_df['LA_Dbulk'] = (LA_D_CPX * m_cpx) + (LA_D_GRT * m_grt) + (LA_D_RUT * m_rut)

    # safe division helper
    def safe_div(a, b):
        a = np.array(a, dtype=float)
        b = np.array(b, dtype=float)
        out = np.full_like(a, np.nan, dtype=float)
        mask = np.isfinite(a) & np.isfinite(b) & (b != 0)
        out[mask] = a[mask] / b[mask]
        return out

    rep_df['MO_F_EM'] = safe_div(rep_df.get('C0_MO', np.nan), rep_df['MO_Dbulk'])
    rep_df['CE_F_EM'] = safe_div(rep_df.get('C0_Ce', np.nan), rep_df['CE_Dbulk'])
    rep_df['W_F_EM'] = safe_div(rep_df.get('C0_W', np.nan), rep_df['W_Dbulk'])
    rep_df['U_F_EM'] = safe_div(rep_df.get('C0_U', np.nan), rep_df['U_Dbulk'])
    rep_df['TH_F_EM'] = safe_div(rep_df.get('C0_Th', np.nan), rep_df['TH_Dbulk'])
    rep_df['NB_F_EM'] = safe_div(rep_df.get('C0_Nb', np.nan), rep_df['NB_Dbulk'])
    rep_df['LA_F_EM'] = safe_div(rep_df.get('C0_La', np.nan), rep_df['LA_Dbulk'])

    # select and return relevant columns
    keep_cols = [
        'rock_id', 'P_bar', 'T_K', 'dFMQ', 'NaCl_wt_pct', 'NaCl_m', 'log10_NaCl',
        'log10_fO2_FMQ', 'log10_fO2_abs',
        'C0_MO', 'C0_Ce', 'C0_W', 'C0_U', 'C0_Th', 'C0_Nb', 'C0_La',
        'MO_D_CPX', 'MO_D_GRT', 'MO_D_RUT', 'W_D_CPX', 'W_D_GRT', 'U_D_CPX', 'U_D_GRT', 'U_D_RUT',
        'MO_Dbulk', 'CE_Dbulk', 'W_Dbulk', 'U_Dbulk', 'TH_Dbulk', 'NB_Dbulk', 'LA_Dbulk',
        'MO_F_EM', 'CE_F_EM', 'W_F_EM', 'U_F_EM', 'TH_F_EM', 'NB_F_EM', 'LA_F_EM'
    ]
    keep = [c for c in keep_cols if c in rep_df.columns]
    out = rep_df[keep].copy()
    out['rock_id'] = out['rock_id'].astype(str)
    return out
