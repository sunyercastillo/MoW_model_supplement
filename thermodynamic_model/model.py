"""
Mo/W–U/Th forward model of Bali et al. (2011, 2012): constants, equations and grid builders.

This is the single source of truth for the supplementary material. The same
definitions drive the Excel workbook (`build_workbook.py`), the DataGraph
tables and the interactive dashboard (`dashboard.py`). Every constant appears
once, in INPUTS; the `Model` class mirrors, line by line, the formulas written
to the workbook.
"""

import os
import math

# ──────────────────────────────────────────────────────────────────────────────
# INPUTS  (name, value-or-formula, unit, note)   — the only place constants live
# ──────────────────────────────────────────────────────────────────────────────
INPUTS = [
    ("Model conditions", None, None, None),
    ("T_C",        1000,                                   "°C",     "Temperature"),
    ("T_K",        "=T_C+273.15",                          "K",      "Temperature in kelvin"),
    ("P_bar",      26100,                                  "bar",    "Pressure (26.1 kbar)"),
    ("NaCl_wt",    10,                                     "wt.%",   "Reference fluid salinity (Rutile and sediment grids)"),
    ("M_NaCl",     58.44,                                  "g/mol",  "Molar mass NaCl"),
    ("NaCl_m",     "=1000*NaCl_wt/(M_NaCl*(100-NaCl_wt))", "mol/kg", "NaCl molality from wt.%"),
    ("log_NaCl",   "=LOG10(NaCl_m)",                       "-",      "log10 NaCl molality"),
    ("Rut_wt",     2.5,                                    "wt.%",   "Reference residual rutile (Salinity and sediment grids)"),
    ("cpx_share",  0.70,                                   "-",      "Clinopyroxene share of the silicate assemblage (cpx:grt = 70:30)"),
    ("x_rut",      "=Rut_wt/100",                          "-",      "Modal rutile"),
    ("x_cpx",      "=cpx_share*(1-x_rut)",                 "-",      "Modal clinopyroxene (0.6825 at 2.5 wt.% rutile)"),
    ("x_grt",      "=(1-cpx_share)*(1-x_rut)",             "-",      "Modal garnet (0.2925 at 2.5 wt.% rutile)"),
    ("f_fluid",    0.01,                                   "-",      "Mass fraction of slab fluid in the metasomatised source (1 %)"),
    ("F_batch",    0,                                      "-",      "Fluid fraction in the batch equilibrium expression C_fluid = C0/(F+(1-F)·D); F = 0 gives C0/D (as in Bali workbook)"),
    ("FMQ buffer (log10 fO2 = A/T + B + C·(P-1)/T)", None, None, None),
    ("FMQ_A",      -25096.3,                               "K",      "Bali et al. (2012) FMQ expression, P in bar"),
    ("FMQ_B",      8.735,                                  "-",      ""),
    ("FMQ_C",      0.11,                                   "K/bar",  ""),
    ("Mo solubility: log10 m_Mo = Mo_a·log fO2 + Mo_b·log m_NaCl + Mo_c·(1000/T) + Mo_d", None, None, None),
    ("Mo_a",       0.44,                                   "-",      "fO2 coefficient. Bali et al. (2012) Appendix C uses 0.435; 0.44 is used here (see audit)"),
    ("Mo_b",       0.42,                                   "-",      "NaCl coefficient"),
    ("Mo_c",       -1.8,                                   "-",      "Temperature coefficient"),
    ("Mo_d",       4.8,                                    "-",      "Intercept"),
    ("M_Mo",       95.95,                                  "g/mol",  "Molar mass Mo"),
    ("Mo_sat_cpx", 40,                                     "µg/g",   "Mo saturation content of cpx (D = c_sat / fluid solubility)"),
    ("Mo_sat_grt", 12,                                     "µg/g",   "Mo saturation content of garnet"),
    ("Mo_sat_rut", 87670,                                  "µg/g",   "Mo saturation content of rutile"),
    ("W solubility: log10 m_W = W_a·log fO2 + W_c·(1000/T) + W_d", None, None, None),
    ("W_a",        0.07,                                   "-",      "fO2 coefficient"),
    ("W_c",        -4.7236,                                "-",      "Temperature coefficient"),
    ("W_d",        4.4271,                                 "-",      "Intercept"),
    ("M_W",        183.84,                                 "g/mol",  "Molar mass W"),
    ("W_sat_cpx",  60,                                     "µg/g",   "W saturation content of cpx"),
    ("W_sat_grt",  12,                                     "µg/g",   "W saturation content of garnet"),
    ("W_D_rut",    1.25,                                   "-",      "Rutile/fluid D for W (constant)"),
    ("U solubility: log10 c_U = U_d + U_a·log fO2 + U_b·m_NaCl  (two intercepts)", None, None, None),
    ("U_a",        0.1433,                                 "-",      "fO2 coefficient (Bali et al. 2011)"),
    ("U_b",        0.594,                                  "kg/mol", "NaCl molality coefficient"),
    ("U_d_sil",    2.681,                                  "-",      "Intercept, silicate hosts (cpx, grt)"),
    ("U_d_rut",    1.7954,                                 "-",      "Intercept, rutile host"),
    ("U_sat_cpx",  11,                                     "µg/g",   "U saturation content of cpx"),
    ("U_sat_grt",  40,                                     "µg/g",   "U saturation content of garnet"),
    ("U_sat_rut",  94,                                     "µg/g",   "U saturation content of rutile"),
    ("Constant mineral/fluid partition coefficients", None, None, None),
    ("Th_D_cpx",   1.19,  "-", "Th cpx/fluid"),
    ("Th_D_grt",   0.61,  "-", "Th grt/fluid"),
    ("Th_D_rut",   0.10,  "-", "Th rut/fluid"),
    ("Nb_D_cpx",   0.172, "-", "Nb cpx/fluid"),
    ("Nb_D_grt",   0.204, "-", "Nb grt/fluid"),
    ("Nb_D_rut",   200,   "-", "Nb rut/fluid"),
    ("La_D_cpx",   1.429, "-", "La cpx/fluid"),
    ("La_D_grt",   0.204, "-", "La grt/fluid"),
    ("La_D_rut",   1.25,  "-", "La rut/fluid"),
    ("Sediment mixture", None, None, None),
    ("carb_frac",  0.6, "-", "Carbonate mass fraction of the mixed sediment (0.6 → 60C:40H)"),
]

# End-member concentrations, µg/g (Table B1). Nb and La are needed only for the
# rutile grid x-axis (Nb/La); their sediment values are not used by any grid.
ELEMENTS = ["Mo", "W", "U", "Th", "Nb", "La"]
ENDMEMBERS = {
    "NMORB": {"Mo": 0.408, "W": 0.038, "U": 0.071, "Th": 0.187, "Nb": 3.507, "La": 3.895},
    "DM":    {"Mo": 0.025, "W": 0.002, "U": 0.005, "Th": 0.014, "Nb": 0.21,  "La": 0.234},
    "Carb":  {"Mo": 1.087, "W": 0.187, "U": 0.684, "Th": 0.222, "Nb": 1.317, "La": 5.254},
    "Hemi":  {"Mo": 5.626, "W": 0.874, "U": 4.117, "Th": 2.656, "Nb": 5.910, "La": 12.862},
}
ENDMEMBER_LABELS = {"NMORB": "N-MORB", "DM": "Depleted mantle", "Carb": "Carbonate", "Hemi": "Hemipelagic"}

DFMQ = [-2, -1, 0, 1, 2, 3]
SED_WT = [0, 1, 2.5, 5, 10]
RUT_WT = [0, 0.25, 0.5, 1, 2.5, 3.5, 5]
NACL_WT = [1, 3, 5, 10, 15, 20]

# Check values of the manuscript grid table at ΔFMQ 0 (computed with Mo_a = 0.435; see docs/audit)
TARGETS_DFMQ0 = {
    # sed wt% : (carb U/Th, carb Mo/W, 6040 U/Th, 6040 Mo/W)
    0:   (0.754334, 8.170971, 0.754334, 8.170971),
    1:   (0.792138, 8.007876, 0.871773, 7.814167),
    2.5: (0.848811, 7.782239, 1.045135, 7.382399),
    5:   (0.943183, 7.449550, 1.326850, 6.849271),
    10:  (1.131609, 6.909933, 1.864670, 6.174343),
}


def param_tag(v):
    """Column-name fragment for a parameter value: 2.5 -> '2_5', 0.25 -> '0_25'."""
    s = ("%g" % v).replace(".", "_")
    return s


# ──────────────────────────────────────────────────────────────────────────────
# Python mirror of the workbook formulas
# ──────────────────────────────────────────────────────────────────────────────
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repository root
MODEL_INPUTS_XLSX = os.path.join(ROOT, "data", "model_inputs.xlsx")
EM_SHEET_COLUMNS = {"N-MORB": "NMORB", "Depleted mantle": "DM", "Carbonate": "Carb", "Hemipelagic": "Hemi"}
GRID_ROWS = {"ΔFMQ (log units)": "dfmq", "Sediment addition (wt.%)": "sed_wt",
             "Residual rutile (wt.%)": "rut_wt", "Fluid salinity (wt.% NaCl)": "nacl_wt"}


def write_model_inputs_template(path=MODEL_INPUTS_XLSX, overwrite=False):
    """
    Write data/model_inputs.xlsx — the user-editable copy of the model inputs
    (Parameters, Endmembers and Grids sheets). Existing files are kept unless
    overwrite=True, so user edits survive a rebuild.
    """
    if os.path.exists(path) and not overwrite:
        return False
    from openpyxl import Workbook
    from openpyxl.styles import Font
    wb = Workbook()
    ws = wb.active
    ws.title = "Parameters"
    ws.append(["Name", "Value", "Unit", "Note"])
    for name, val, unit, note in INPUTS:
        if val is None:
            ws.append([name, None, None, None])
            ws.cell(row=ws.max_row, column=1).font = Font(bold=True)
        elif isinstance(val, str):
            continue  # derived quantities are recomputed, not edited
        else:
            ws.append([name, val, unit, note])
    ws.column_dimensions["A"].width = 14
    ws.column_dimensions["D"].width = 90
    ws2 = wb.create_sheet("Endmembers")
    ws2.append(["Element"] + list(EM_SHEET_COLUMNS))
    for e in ELEMENTS:
        ws2.append([e] + [ENDMEMBERS[k][e] for k in EM_SHEET_COLUMNS.values()])
    ws2.append([])
    ws2.append(["Concentrations in µg/g. The mixed sediment is carb_frac·Carbonate + (1 − carb_frac)·Hemipelagic "
                "(carb_frac in Parameters); the dashboard slider overrides carb_frac."])
    ws3 = wb.create_sheet("Grids")
    ws3.append(["Axis", "Values →"])
    for label, attr in GRID_ROWS.items():
        ws3.append([label] + {"dfmq": DFMQ, "sed_wt": SED_WT, "rut_wt": RUT_WT, "nacl_wt": NACL_WT}[attr])
    ws3.append([])
    ws3.append(["Values are read left to right; leave cells blank to shorten a list. Sediment list must contain 0."])
    for w in (ws, ws2, ws3):
        for c in w[1]:
            c.font = Font(bold=True)
    wb.save(path)
    return True


def load_model_inputs(path=MODEL_INPUTS_XLSX):
    """Return (constants, endmembers, grids) read from data/model_inputs.xlsx, or Python defaults if absent."""
    c = {name: float(val) for name, val, _, _ in INPUTS if val is not None and not isinstance(val, str)}
    em = {k: dict(v) for k, v in ENDMEMBERS.items()}
    grids = {"dfmq": list(DFMQ), "sed_wt": list(SED_WT), "rut_wt": list(RUT_WT), "nacl_wt": list(NACL_WT)}
    if not os.path.exists(path):
        return c, em, grids, None
    from openpyxl import load_workbook
    wb = load_workbook(path, data_only=True)
    ws = wb["Parameters"]
    for row in ws.iter_rows(min_row=2, values_only=True):
        if row[0] and row[1] is not None and str(row[0]) in c:
            c[str(row[0])] = float(row[1])
    ws = wb["Endmembers"]
    header = [cell.value for cell in ws[1]]
    for row in ws.iter_rows(min_row=2, values_only=True):
        if row[0] in ELEMENTS:
            for j, h in enumerate(header[1:], start=1):
                if h in EM_SHEET_COLUMNS and row[j] is not None:
                    em[EM_SHEET_COLUMNS[h]][row[0]] = float(row[j])
    ws = wb["Grids"]
    for row in ws.iter_rows(min_row=2, values_only=True):
        if row[0] in GRID_ROWS:
            vals = [float(v) for v in row[1:] if v is not None and v != ""]
            vals = [int(v) if v == int(v) else v for v in vals]
            if vals:
                grids[GRID_ROWS[row[0]]] = vals
    return c, em, grids, path


class Model:
    """
    Forward model at fixed conditions. Constants, end-members and grid axes are
    taken from data/model_inputs.xlsx when that file exists (user-editable), and
    from the Python defaults in this module otherwise.
    """

    def __init__(self, mo_a=None, inputs_path=MODEL_INPUTS_XLSX):
        self.c, self.em, grids, self.inputs_file = load_model_inputs(inputs_path)
        self.dfmq, self.sed_wt, self.rut_wt, self.nacl_wt = (grids["dfmq"], grids["sed_wt"], grids["rut_wt"], grids["nacl_wt"])
        if mo_a is not None:
            self.c["Mo_a"] = mo_a
        c = self.c
        c["T_K"] = c["T_C"] + 273.15
        c["NaCl_m"] = 1000 * c["NaCl_wt"] / (c["M_NaCl"] * (100 - c["NaCl_wt"]))
        c["log_NaCl"] = math.log10(c["NaCl_m"])
        c["x_rut"] = c["Rut_wt"] / 100
        c["x_cpx"] = c["cpx_share"] * (1 - c["x_rut"])
        c["x_grt"] = (1 - c["cpx_share"]) * (1 - c["x_rut"])
        self.mix = self.mixed_sediment(c["carb_frac"])

    def row(self, dfmq, protolith, nacl_wt=None, rut_wt=None):
        c = self.c
        nacl_wt = c["NaCl_wt"] if nacl_wt is None else nacl_wt
        rut_wt = c["Rut_wt"] if rut_wt is None else rut_wt
        r = {}
        r["NaCl_m"] = 1000 * nacl_wt / (c["M_NaCl"] * (100 - nacl_wt))
        r["log_NaCl"] = math.log10(r["NaCl_m"])
        r["x_rut"] = rut_wt / 100
        r["x_cpx"] = c["cpx_share"] * (1 - r["x_rut"])
        r["x_grt"] = (1 - c["cpx_share"]) * (1 - r["x_rut"])
        r["logfO2_FMQ"] = c["FMQ_A"] / c["T_K"] + c["FMQ_B"] + c["FMQ_C"] * ((c["P_bar"] - 1) / c["T_K"])
        r["logfO2"] = r["logfO2_FMQ"] + dfmq
        r["Mo_sol_m"] = 10 ** (c["Mo_a"] * r["logfO2"] + c["Mo_b"] * r["log_NaCl"] + c["Mo_c"] * (1000 / c["T_K"]) + c["Mo_d"])
        r["Mo_sol_ugg"] = r["Mo_sol_m"] * c["M_Mo"] * 1000
        r["W_sol_m"] = 10 ** (c["W_a"] * r["logfO2"] + c["W_c"] * (1000 / c["T_K"]) + c["W_d"])
        r["W_sol_ugg"] = r["W_sol_m"] * c["M_W"] * 1000
        r["U_sol_sil"] = 10 ** (c["U_d_sil"] + c["U_a"] * r["logfO2"] + c["U_b"] * r["NaCl_m"])
        r["U_sol_rut"] = 10 ** (c["U_d_rut"] + c["U_a"] * r["logfO2"] + c["U_b"] * r["NaCl_m"])
        x_cpx, x_grt, x_rut = r["x_cpx"], r["x_grt"], r["x_rut"]
        r["D_Mo_cpx"] = c["Mo_sat_cpx"] / r["Mo_sol_ugg"]
        r["D_Mo_grt"] = c["Mo_sat_grt"] / r["Mo_sol_ugg"]
        r["D_Mo_rut"] = c["Mo_sat_rut"] / r["Mo_sol_ugg"]
        r["D_W_cpx"] = c["W_sat_cpx"] / r["W_sol_ugg"]
        r["D_W_grt"] = c["W_sat_grt"] / r["W_sol_ugg"]
        r["D_U_cpx"] = c["U_sat_cpx"] / r["U_sol_sil"]
        r["D_U_grt"] = c["U_sat_grt"] / r["U_sol_sil"]
        r["D_U_rut"] = c["U_sat_rut"] / r["U_sol_rut"]
        D = {
            "Mo": r["D_Mo_cpx"] * x_cpx + r["D_Mo_grt"] * x_grt + r["D_Mo_rut"] * x_rut,
            "W":  r["D_W_cpx"] * x_cpx + r["D_W_grt"] * x_grt + c["W_D_rut"] * x_rut,
            "U":  r["D_U_cpx"] * x_cpx + r["D_U_grt"] * x_grt + r["D_U_rut"] * x_rut,
            "Th": c["Th_D_cpx"] * x_cpx + c["Th_D_grt"] * x_grt + c["Th_D_rut"] * x_rut,
            "Nb": c["Nb_D_cpx"] * x_cpx + c["Nb_D_grt"] * x_grt + c["Nb_D_rut"] * x_rut,
            "La": c["La_D_cpx"] * x_cpx + c["La_D_grt"] * x_grt + c["La_D_rut"] * x_rut,
        }
        F, f = c["F_batch"], c["f_fluid"]
        for e in ELEMENTS:
            r["D_" + e] = D[e]
            r["Prot_" + e] = protolith[e]
            r["Fluid_" + e] = protolith[e] / (F + (1 - F) * D[e])
            r["Src_" + e] = f * r["Fluid_" + e] + (1 - f) * self.em["DM"][e]
        r["UTh"] = r["Src_U"] / r["Src_Th"]
        r["MoW"] = r["Src_Mo"] / r["Src_W"]
        r["NbLa"] = r["Src_Nb"] / r["Src_La"]
        r["Fluid_NbLa"] = r["Fluid_Nb"] / r["Fluid_La"]
        r["Fluid_MoW"] = r["Fluid_Mo"] / r["Fluid_W"]
        return r

    def mixed_sediment(self, carb_frac):
        """Sediment end-member with carbonate mass fraction carb_frac (0.6 → 60C:40H)."""
        return {e: carb_frac * self.em["Carb"][e] + (1 - carb_frac) * self.em["Hemi"][e] for e in ELEMENTS}

    def protolith(self, s, sediment):
        return {e: (1 - s) * self.em["NMORB"][e] + s * sediment[e] for e in ELEMENTS}

    def sediment_grid(self, sediment):
        return [(d, w, self.row(d, self.protolith(w / 100, sediment))) for d in self.dfmq for w in self.sed_wt]

    def rutile_grid(self):
        return [(d, w, self.row(d, self.em["NMORB"], rut_wt=w)) for d in self.dfmq for w in self.rut_wt]

    def salinity_grid(self):
        return [(d, w, self.row(d, self.em["NMORB"], nacl_wt=w)) for d in self.dfmq for w in self.nacl_wt]


def dg_table(grid_rows, xkey, iso_vals, param_vals, ykey="MoW"):
    """Python-mirror version of the DataGraph table: list of column (header, values)."""
    n = len(param_vals)
    n_dfmq = len(grid_rows) // n
    X = [r[xkey] for _, _, r in grid_rows]
    Y = [r[ykey] for _, _, r in grid_rows]
    idx = lambda i_d, i_p: i_d * n + i_p
    cols = [("FMQ", [d for d, _, _ in grid_rows]), ("PARAM", [w for _, w, _ in grid_rows]),
            ("X", X), ("MoW", Y), ("Band_mask", [1] * (2 * n + 1))]
    for b in range(n_dfmq - 1):
        seq = [(b, i) for i in range(n)] + [(b + 1, i) for i in reversed(range(n))] + [(b, 0)]
        cols.append((f"Band{b + 1}_X", [X[idx(*s)] for s in seq]))
        cols.append((f"Band{b + 1}_Y", [Y[idx(*s)] for s in seq]))
    for v in iso_vals:
        i_p = param_vals.index(v)
        cols.append((f"ISO{param_tag(v)}_X", [X[idx(i_d, i_p)] for i_d in range(n_dfmq)]))
        cols.append((f"ISO{param_tag(v)}_Y", [Y[idx(i_d, i_p)] for i_d in range(n_dfmq)]))
    return cols


def fmt_num(x):
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        if x == int(x) and abs(x) < 1e6:
            return str(int(x))
        return "%.10g" % x
    return str(x)


def tsv_text(headers, cols):
    """Tab-separated text of a DataGraph table (columns of unequal length, blanks after the end)."""
    nrow = max(len(v) for _, v in cols)
    lines = ["\t".join(headers)]
    for i in range(nrow):
        lines.append("\t".join(fmt_num(v[i]) if i < len(v) else "" for _, v in cols))
    return "\n".join(lines) + "\n"


def grid_dataframe_rows(grid_rows, param_label):
    """Flat list of dicts (one per grid point) with every intermediate, for display tables."""
    out = []
    for dfmq, w, r in grid_rows:
        d = {"ΔFMQ": dfmq, param_label: w}
        d.update(r)
        out.append(d)
    return out
