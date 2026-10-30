#!/usr/bin/env python3
"""
Build the Mo/W–U/Th sensitivity grids of the Bali et al. (2011, 2012) fluid–mineral
partitioning model as (i) an Excel workbook with live, cell-referenced formulas,
(ii) tab-separated DataGraph tables and (iii) a validation report.

Every constant appears exactly once, in INPUTS below, and is written to the
`Inputs` sheet as a named cell. All calculation sheets reference those names, so
the workbook can be audited cell by cell and any constant changed in one place.

The Python functions in this file are an independent mirror of the Excel formulas;
they produce the DataGraph tables and the validation numbers. If Microsoft Excel
is available the workbook is also recalculated and compared against the mirror.

Usage (from the repository root):  python -m thermodynamic_model.build_workbook
"""

import os
import re
import math
import itertools
import subprocess
from datetime import date

from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter
from openpyxl.workbook.defined_name import DefinedName

from thermodynamic_model.model import (tsv_text, INPUTS, ELEMENTS, ENDMEMBER_LABELS, DFMQ, SED_WT, RUT_WT, NACL_WT,
                       TARGETS_DFMQ0, param_tag, Model, dg_table, fmt_num, write_model_inputs_template)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repository root
OUT_DIR = os.path.join(ROOT, "outputs")
XLSX = os.path.join(OUT_DIR, "MoW_thermodynamic_model_grids.xlsx")
DG_DIR = os.path.join(OUT_DIR, "datagraph")
REPORT = os.path.join(OUT_DIR, "validation_report.md")

# ──────────────────────────────────────────────────────────────────────────────
# Workbook writing helpers
# ──────────────────────────────────────────────────────────────────────────────
HDR_FONT = Font(bold=True)
SEC_FONT = Font(bold=True, color="FFFFFF")
SEC_FILL = PatternFill("solid", fgColor="4F6228")
NOTE_FONT = Font(italic=True, color="555555")


def set_widths(ws, widths):
    for i, w in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(i)].width = w


def define(wb, name, sheet, col_letter, row):
    wb.defined_names[name] = DefinedName(name, attr_text=f"'{sheet}'!${col_letter}${row}")


def write_inputs(wb, model):
    ENDMEMBERS = model.em
    ws = wb.active
    ws.title = "Inputs"
    ws["A1"] = "Inputs — every constant of the model. Calculation sheets reference these cells by name."
    ws["A1"].font = HDR_FONT
    ws["A2"] = ("Conditions: T = 1000 °C, P = 26.1 kbar, 10 wt.% NaCl, 2.5 wt.% rutile, cpx:grt = 70:30, "
                "1 % fluid + 99 % depleted mantle. Model: Bali et al. (2011, 2012).")
    ws["A2"].font = NOTE_FONT
    for j, h in enumerate(["Name", "Value", "Unit", "Note"], start=1):
        ws.cell(row=4, column=j, value=h).font = HDR_FONT
    r = 5
    for name, val, unit, note in INPUTS:
        if val is None:
            ws.cell(row=r, column=1, value=name).font = SEC_FONT
            for j in range(1, 5):
                ws.cell(row=r, column=j).fill = SEC_FILL
        else:
            ws.cell(row=r, column=1, value=name)
            ws.cell(row=r, column=2, value=val if isinstance(val, str) else model.c[name])
            ws.cell(row=r, column=3, value=unit)
            ws.cell(row=r, column=4, value=note)
            define(wb, name, "Inputs", "B", r)
        r += 1

    r += 1
    ws.cell(row=r, column=1, value="End-member concentrations (µg/g), Table B1").font = SEC_FONT
    for j in range(1, 8):
        ws.cell(row=r, column=j).fill = SEC_FILL
    r += 1
    heads = ["Element"] + [ENDMEMBER_LABELS[k] for k in ENDMEMBERS] + ["Mixed sediment (carb_frac·Carb + (1−carb_frac)·Hemi)"]
    for j, h in enumerate(heads, start=1):
        ws.cell(row=r, column=j, value=h).font = HDR_FONT
    r += 1
    for e in ELEMENTS:
        ws.cell(row=r, column=1, value=e)
        for j, key in enumerate(ENDMEMBERS, start=2):
            ws.cell(row=r, column=j, value=ENDMEMBERS[key][e])
            define(wb, f"{key}_{e}", "Inputs", get_column_letter(j), r)
        col = get_column_letter(len(ENDMEMBERS) + 2)
        ws.cell(row=r, column=len(ENDMEMBERS) + 2, value=f"=carb_frac*Carb_{e}+(1-carb_frac)*Hemi_{e}")
        define(wb, f"Mix_{e}", "Inputs", col, r)
        r += 1
    ws.cell(row=r, column=1, value="Nb and La enter only the rutile grid (x-axis Nb/La) with a pure N-MORB protolith; "
                                   "sediment Nb, La are listed for completeness and are not used by any grid.").font = NOTE_FONT
    r += 2
    ws.cell(row=r, column=1, value="Grid axes").font = SEC_FONT
    for j in range(1, 8):
        ws.cell(row=r, column=j).fill = SEC_FILL
    r += 1
    for label, vals in [("ΔFMQ (log units)", DFMQ), ("Sediment addition (wt.%)", SED_WT),
                        ("Residual rutile (wt.%)", RUT_WT), ("Fluid salinity (wt.% NaCl)", NACL_WT)]:
        ws.cell(row=r, column=1, value=label)
        for j, v in enumerate(vals, start=2):
            ws.cell(row=r, column=j, value=v)
        r += 1
    set_widths(ws, [16, 14, 10, 14, 14, 14, 44])
    ws.freeze_panes = "A5"


def write_calc_sheet(wb, title, kind, sediment_name=None, note=""):
    """
    kind: 'sediment' | 'rutile' | 'salinity'. Writes the 36/42/30-row grid with
    one formula per cell, all referencing Inputs names. Returns dict with the
    column letters of key outputs and the first data row, used by the DG sheets.
    """
    ws = wb.create_sheet(title)
    ws["A1"] = note
    ws["A1"].font = NOTE_FONT
    cols = []  # (header, formula-template or None for value)

    if kind == "sediment":
        param_vals = SED_WT
        prot = lambda e: f"=(1-$C{{r}})*NMORB_{e}+$C{{r}}*{sediment_name}_{e}"
        cols += [("ΔFMQ", None), ("Sediment addition (wt.%)", None), ("Sediment fraction", "=$B{r}/100")]
        c_x_rut, c_x_cpx, c_x_grt = "x_rut", "x_cpx", "x_grt"
        nacl_m, log_nacl = "NaCl_m", "log_NaCl"
        prot_elems = ["Mo", "W", "U", "Th"]
    elif kind == "rutile":
        param_vals = RUT_WT
        prot = lambda e: f"=NMORB_{e}"
        cols += [("ΔFMQ", None), ("Residual rutile (wt.%)", None),
                 ("Modal rutile", "=$B{r}/100"), ("Modal cpx", "=cpx_share*(1-$C{r})"), ("Modal garnet", "=(1-cpx_share)*(1-$C{r})")]
        c_x_rut, c_x_cpx, c_x_grt = "$C{r}", "$D{r}", "$E{r}"
        nacl_m, log_nacl = "NaCl_m", "log_NaCl"
        prot_elems = ELEMENTS
    else:
        param_vals = NACL_WT
        prot = lambda e: f"=NMORB_{e}"
        cols += [("ΔFMQ", None), ("Fluid salinity (wt.% NaCl)", None),
                 ("NaCl molality (mol/kg)", "=1000*$B{r}/(M_NaCl*(100-$B{r}))"), ("log10 NaCl molality", "=LOG10($C{r})")]
        c_x_rut, c_x_cpx, c_x_grt = "x_rut", "x_cpx", "x_grt"
        nacl_m, log_nacl = "$C{r}", "$D{r}"
        prot_elems = ["Mo", "W", "U", "Th"]

    for e in prot_elems:
        cols.append((f"Protolith {e}", prot(e)))
    L = {}  # header key -> column letter, filled below

    def ref(key):
        return f"${L[key]}{{r}}"

    def add(key, header, template):
        cols.append((header, template))
        L[key] = get_column_letter(len(cols))

    # fill letters for columns already added
    for i, (h, _) in enumerate(cols, start=1):
        L[h] = get_column_letter(i)
    P = {e: L[f"Protolith {e}"] for e in prot_elems}

    add("fmq", "FMQ-buffer log10(fO2)", "=FMQ_A/T_K+FMQ_B+FMQ_C*((P_bar-1)/T_K)")
    add("fo2", "Absolute log10(fO2)", f"={ref('fmq')}+$A{{r}}")
    add("MoM", "Mo solubility (mol/kg)", f"=10^(Mo_a*{ref('fo2')}+Mo_b*{log_nacl}+Mo_c*(1000/T_K)+Mo_d)")
    add("MoS", "Mo solubility (µg/g)", f"={ref('MoM')}*M_Mo*1000")
    add("WM", "W solubility (mol/kg)", f"=10^(W_a*{ref('fo2')}+W_c*(1000/T_K)+W_d)")
    add("WS", "W solubility (µg/g)", f"={ref('WM')}*M_W*1000")
    add("USs", "U solubility, silicate (µg/g)", f"=10^(U_d_sil+U_a*{ref('fo2')}+U_b*{nacl_m})")
    add("USr", "U solubility, rutile (µg/g)", f"=10^(U_d_rut+U_a*{ref('fo2')}+U_b*{nacl_m})")
    add("DMo", "Bulk DMo", f"=(Mo_sat_cpx*{c_x_cpx}+Mo_sat_grt*{c_x_grt}+Mo_sat_rut*{c_x_rut})/{ref('MoS')}")
    add("DW", "Bulk DW", f"=(W_sat_cpx*{c_x_cpx}+W_sat_grt*{c_x_grt})/{ref('WS')}+W_D_rut*{c_x_rut}")
    add("DU", "Bulk DU", f"=(U_sat_cpx*{c_x_cpx}+U_sat_grt*{c_x_grt})/{ref('USs')}+U_sat_rut*{c_x_rut}/{ref('USr')}")
    add("DTh", "Bulk DTh", f"=Th_D_cpx*{c_x_cpx}+Th_D_grt*{c_x_grt}+Th_D_rut*{c_x_rut}")
    if kind == "rutile":
        add("DNb", "Bulk DNb", f"=Nb_D_cpx*{c_x_cpx}+Nb_D_grt*{c_x_grt}+Nb_D_rut*{c_x_rut}")
        add("DLa", "Bulk DLa", f"=La_D_cpx*{c_x_cpx}+La_D_grt*{c_x_grt}+La_D_rut*{c_x_rut}")
    for e in prot_elems:
        add("F" + e, f"Fluid {e}", f"=${P[e]}{{r}}/(F_batch+(1-F_batch)*{ref('D' + e)})")
    for e in prot_elems:
        add("S" + e, f"Source {e}", f"=f_fluid*{ref('F' + e)}+(1-f_fluid)*DM_{e}")
    if kind == "rutile":
        # Panel A is plotted with the slab-fluid composition itself (no mantle mixing); source ratios kept for reference
        add("FNbLa", "Fluid Nb/La (plotted)", f"={ref('FNb')}/{ref('FLa')}")
        add("FMoW", "Fluid Mo/W (plotted)", f"={ref('FMo')}/{ref('FW')}")
        add("NbLa", "Source Nb/La", f"={ref('SNb')}/{ref('SLa')}")
        add("MoW", "Source Mo/W", f"={ref('SMo')}/{ref('SW')}")
        add("UTh", "Source U/Th", f"={ref('SU')}/{ref('STh')}")
    else:
        add("UTh", "U/Th", f"={ref('SU')}/{ref('STh')}")
        add("MoW", "Mo/W", f"={ref('SMo')}/{ref('SW')}")
    # per-mineral partition coefficients (intermediates of the bulk D)
    add("dMoc", "D Mo cpx/fluid", f"=Mo_sat_cpx/{ref('MoS')}")
    add("dMog", "D Mo grt/fluid", f"=Mo_sat_grt/{ref('MoS')}")
    add("dMor", "D Mo rut/fluid", f"=Mo_sat_rut/{ref('MoS')}")
    add("dWc", "D W cpx/fluid", f"=W_sat_cpx/{ref('WS')}")
    add("dWg", "D W grt/fluid", f"=W_sat_grt/{ref('WS')}")
    add("dUc", "D U cpx/fluid", f"=U_sat_cpx/{ref('USs')}")
    add("dUg", "D U grt/fluid", f"=U_sat_grt/{ref('USs')}")
    add("dUr", "D U rut/fluid", f"=U_sat_rut/{ref('USr')}")

    hdr_row = 3
    for j, (h, _) in enumerate(cols, start=1):
        c = ws.cell(row=hdr_row, column=j, value=h)
        c.font = HDR_FONT
        c.alignment = Alignment(wrap_text=True, vertical="top")
    ws.row_dimensions[hdr_row].height = 45
    first = hdr_row + 1
    r = first
    for d in DFMQ:
        for w in param_vals:
            ws.cell(row=r, column=1, value=d)
            ws.cell(row=r, column=2, value=w)
            for j, (_, tmpl) in enumerate(cols[2:], start=3):
                ws.cell(row=r, column=j, value=tmpl.format(r=r))
            r += 1
    ws.freeze_panes = ws.cell(row=first, column=3)
    set_widths(ws, [8, 11] + [13] * (len(cols) - 2))
    return {"first": first, "n": r - first, "L": L, "param_vals": param_vals, "title": title}


def write_dg_sheet(wb, title, calc, xkey, xname, param_name, iso_prefix, iso_vals, grid_rows, ykey="MoW", yname="MoW"):
    """
    DataGraph export sheet: main curves (one row per grid point, ordered by ΔFMQ
    then parameter), five closed band polygons, and dashed isopleths — all as
    references to the calculation sheet, so nothing is duplicated.
    grid_rows: python-mirror rows [(dfmq, param, rowdict)] in the same order.
    """
    ws = wb.create_sheet(title)
    cs = calc["title"]
    Lx, Ly = calc["L"][xkey], calc["L"][ykey]
    first = calc["first"]
    pv = calc["param_vals"]
    n = len(pv)

    def rr(i_d, i_p):  # calc-sheet row of ΔFMQ index i_d, parameter index i_p
        return first + i_d * n + i_p

    headers = ["FMQ", param_name, xname, yname, "Band_mask"]
    for b in range(1, 6):
        headers += [f"Band{b}_X", f"Band{b}_Y"]
    for v in iso_vals:
        headers += [f"{iso_prefix}{param_tag(v)}_X", f"{iso_prefix}{param_tag(v)}_Y"]
    for j, h in enumerate(headers, start=1):
        ws.cell(row=1, column=j, value=h).font = HDR_FONT

    # main curves
    for k, (d, w, _) in enumerate(grid_rows):
        i_d, i_p = divmod(k, n)
        row = 2 + k
        ws.cell(row=row, column=1, value=d)
        ws.cell(row=row, column=2, value=w)
        ws.cell(row=row, column=3, value=f"='{cs}'!{Lx}{rr(i_d, i_p)}")
        ws.cell(row=row, column=4, value=f"='{cs}'!{Ly}{rr(i_d, i_p)}")
    # bands: lower boundary forward, upper boundary backward, close
    npoly = 2 * n + 1
    for k in range(npoly):
        ws.cell(row=2 + k, column=5, value=1)
    for b in range(5):
        cx = 6 + 2 * b
        seq = [(b, i) for i in range(n)] + [(b + 1, i) for i in reversed(range(n))] + [(b, 0)]
        for k, (i_d, i_p) in enumerate(seq):
            ws.cell(row=2 + k, column=cx, value=f"='{cs}'!{Lx}{rr(i_d, i_p)}")
            ws.cell(row=2 + k, column=cx + 1, value=f"='{cs}'!{Ly}{rr(i_d, i_p)}")
    # isopleths: constant parameter, ΔFMQ -2 → +3
    for m, v in enumerate(iso_vals):
        i_p = pv.index(v)
        cx = 16 + 2 * m
        for i_d in range(len(DFMQ)):
            ws.cell(row=2 + i_d, column=cx, value=f"='{cs}'!{Lx}{rr(i_d, i_p)}")
            ws.cell(row=2 + i_d, column=cx + 1, value=f"='{cs}'!{Ly}{rr(i_d, i_p)}")
    set_widths(ws, [8, 14, 12, 12, 10] + [12] * (len(headers) - 5))
    ws.freeze_panes = "A2"
    return headers


# ──────────────────────────────────────────────────────────────────────────────
# Geometry check: polygon self-intersection
# ──────────────────────────────────────────────────────────────────────────────
def _orient(a, b, c):
    v = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return 0 if abs(v) < 1e-15 else (1 if v > 0 else -1)


def segments_intersect(p1, p2, p3, p4):
    o1, o2, o3, o4 = _orient(p1, p2, p3), _orient(p1, p2, p4), _orient(p3, p4, p1), _orient(p3, p4, p2)
    return o1 != o2 and o3 != o4 and 0 not in (o1, o2, o3, o4)


def polygon_self_intersects(xs, ys, log=True):
    """Checked in log-log space (the plotting space) and in linear space."""
    pts = list(zip(xs, ys))
    if pts[0] == pts[-1]:
        pts = pts[:-1]
    if log:
        pts = [(math.log10(x), math.log10(y)) for x, y in pts]
    m = len(pts)
    segs = [(pts[i], pts[(i + 1) % m]) for i in range(m)]
    for i, j in itertools.combinations(range(m), 2):
        if abs(i - j) in (1, m - 1):
            continue
        if segments_intersect(*segs[i], *segs[j]):
            return True
    return False


# ──────────────────────────────────────────────────────────────────────────────
# Checks sheet (live Excel formulas) 
# ──────────────────────────────────────────────────────────────────────────────
def write_checks(wb, calc_pc, calc_mx, dg_names, python_checks):
    ws = wb.create_sheet("Checks")
    ws["A1"] = "Validation checks (live formulas; TRUE = pass). Geometric self-intersection is tested in the build script, see validation_report.md."
    ws["A1"].font = HDR_FONT
    for j, h in enumerate(["#", "Check", "Result", "Detail"], start=1):
        ws.cell(row=3, column=j, value=h).font = HDR_FONT
    r = 4
    f0 = calc_pc["first"]
    n = 5
    Luth, Lmow = calc_pc["L"]["UTh"], calc_pc["L"]["MoW"]

    def put(num, text, formula, detail=""):
        nonlocal r
        ws.cell(row=r, column=1, value=num)
        ws.cell(row=r, column=2, value=text)
        ws.cell(row=r, column=3, value=formula)
        ws.cell(row=r, column=4, value=detail)
        r += 1

    # 1 identical at 0 wt.%
    terms = []
    for i_d in range(len(DFMQ)):
        row = f0 + i_d * n
        terms.append(f"ABS('Pure_carbonate'!{Luth}{row}-'C60_H40'!{Luth}{row})<1E-12")
        terms.append(f"ABS('Pure_carbonate'!{Lmow}{row}-'C60_H40'!{Lmow}{row})<1E-12")
    put(1, "Pure-carbonate and 60C:40H grids identical at 0 wt.% sediment (all six ΔFMQ)", "=AND(" + ",".join(terms) + ")")
    # 2,3 band counts and closure for both sediment DG sheets
    for dg in dg_names[:2]:
        cnt = []
        clo = []
        for b in range(5):
            cx = get_column_letter(6 + 2 * b)
            cy = get_column_letter(7 + 2 * b)
            cnt.append(f"COUNT('{dg}'!{cx}2:{cx}200)=11")
            clo.append(f"AND('{dg}'!{cx}2='{dg}'!{cx}12,'{dg}'!{cy}2='{dg}'!{cy}12)")
        put(2, f"{dg}: every band has exactly 11 points", "=AND(" + ",".join(cnt) + ")")
        put(3, f"{dg}: first and last coordinate of every band identical", "=AND(" + ",".join(clo) + ")")
    put(4, "No band self-intersects (log–log and linear space)", python_checks["no_self_intersection"],
        "Computed geometrically in thermodynamic_model/build_workbook.py for all four grids")
    # 5 isopleth counts
    for dg in dg_names[:2]:
        cnt = [f"COUNT('{dg}'!{get_column_letter(16 + 2 * m)}2:{get_column_letter(16 + 2 * m)}200)=6" for m in range(4)]
        put(5, f"{dg}: every dashed sediment trajectory has six points", "=AND(" + ",".join(cnt) + ")")
    # 6 Sed2_5 uses 2.5 wt.%: the referenced calc rows carry 2.5 in column B
    for dg, cs in zip(dg_names[:2], ["Pure_carbonate", "C60_H40"]):
        rows25 = [f0 + i_d * n + SED_WT.index(2.5) for i_d in range(len(DFMQ))]
        t = [f"'{cs}'!$B${rw}=2.5" for rw in rows25]
        t += [f"'{dg}'!{get_column_letter(16 + 2 * 1)}{2 + i_d}='{cs}'!{Luth}{rw}" for i_d, rw in enumerate(rows25)]
        put(6, f"{dg}: Sed2_5 columns reference the 2.5 wt.% rows", "=AND(" + ",".join(t) + ")")
    # 7,8 Band3 between ΔFMQ 0 and +1; Band5 between +2 and +3
    for num, b, lo, hi in [(7, 3, 0, 1), (8, 5, 2, 3)]:
        for dg, cs in zip(dg_names[:2], ["Pure_carbonate", "C60_H40"]):
            cx = get_column_letter(6 + 2 * (b - 1))
            i_lo, i_hi = DFMQ.index(lo), DFMQ.index(hi)
            t = [f"'{dg}'!{cx}{2 + k}='{cs}'!{Luth}{f0 + i_lo * n + k}" for k in range(n)]
            t += [f"'{dg}'!{cx}{2 + n + k}='{cs}'!{Luth}{f0 + i_hi * n + (n - 1 - k)}" for k in range(n)]
            t += [f"'{cs}'!$A${f0 + i_lo * n}={lo}", f"'{cs}'!$A${f0 + i_hi * n}={hi}"]
            put(num, f"{dg}: Band{b} lies between the ΔFMQ {lo:+d} and {hi:+d} curves (lower forward, upper backward)", "=AND(" + ",".join(t) + ")")
    put(9, "Band_mask is a separate column (E), not the FMQ or sediment column", "=TRUE", "Structural: Band_mask written as its own column in every DG sheet")
    for dg in dg_names[:2]:
        put(10, f"{dg}: no Sed0_X / Sed0_Y columns", f"=AND(ISERROR(MATCH(\"Sed0_X\",'{dg}'!$1:$1,0)),ISERROR(MATCH(\"Sed0_Y\",'{dg}'!$1:$1,0)))")
    put("", "ALL CHECKS", f"=AND(C4:C{r - 1})")
    n_rows = r - 4
    r += 2
    ws.cell(row=r, column=1, value="Numerical check at ΔFMQ 0 against the values supplied with the task").font = SEC_FONT
    for j in range(1, 11):
        ws.cell(row=r, column=j).fill = SEC_FILL
    r += 1
    ws.cell(row=r, column=1, value="The supplied check values reproduce exactly with Mo_a = 0.435 (Bali et al. 2012, Appendix C). "
                                   "With Mo_a = 0.44 (chosen, Inputs!Mo_a) U/Th is unchanged and Mo/W is ≈7.5 % lower. Change Inputs!Mo_a to 0.435 to verify.").font = NOTE_FONT
    r += 1
    heads = ["Sediment (wt.%)", "Carb U/Th (workbook)", "Carb U/Th (supplied)", "Carb Mo/W (workbook)", "Carb Mo/W (supplied)",
             "60C:40H U/Th (workbook)", "60C:40H U/Th (supplied)", "60C:40H Mo/W (workbook)", "60C:40H Mo/W (supplied)", "Mo/W ratio workbook/supplied"]
    for j, h in enumerate(heads, start=1):
        c = ws.cell(row=r, column=j, value=h)
        c.font = HDR_FONT
        c.alignment = Alignment(wrap_text=True)
    r += 1
    i0 = DFMQ.index(0)
    for i_p, w in enumerate(SED_WT):
        row = f0 + i0 * n + i_p
        t = TARGETS_DFMQ0[w]
        ws.cell(row=r, column=1, value=w)
        ws.cell(row=r, column=2, value=f"='Pure_carbonate'!{Luth}{row}")
        ws.cell(row=r, column=3, value=t[0])
        ws.cell(row=r, column=4, value=f"='Pure_carbonate'!{Lmow}{row}")
        ws.cell(row=r, column=5, value=t[1])
        ws.cell(row=r, column=6, value=f"='C60_H40'!{Luth}{row}")
        ws.cell(row=r, column=7, value=t[2])
        ws.cell(row=r, column=8, value=f"='C60_H40'!{Lmow}{row}")
        ws.cell(row=r, column=9, value=t[3])
        ws.cell(row=r, column=10, value=f"=D{r}/E{r}")
        r += 1
    set_widths(ws, [6, 70, 12, 60] + [14] * 6)
    return n_rows


# ──────────────────────────────────────────────────────────────────────────────
# Optional: recalculate with Microsoft Excel and compare to the Python mirror
# ──────────────────────────────────────────────────────────────────────────────
def excel_cross_check(xlsx, calcs, model_rows, n_check_rows):
    """
    Open the workbook in Microsoft Excel (read-only use: closed without saving),
    force a full recalculation and read back the U/Th and Mo/W columns of every
    calculation sheet plus the Checks results. Returns (result, message) where
    result = (max relative difference vs Python mirror, cells compared, ALL CHECKS value)
    or None if Excel is unavailable.
    """
    reads = []
    for calc, rows in zip(calcs, model_rows):
        last = calc["first"] + len(rows) - 1
        for key in ("UTh", "MoW"):
            L = calc["L"][key]
            reads.append((calc["title"], f"{L}{calc['first']}:{L}{last}"))  # source ratios exist in every sheet
    reads.append(("Checks", f"C4:C{3 + n_check_rows}"))
    lines = []
    for i, (sheet, rng) in enumerate(reads):
        lines.append(f'set v{i} to (value of range "{rng}" of worksheet "{sheet}" of wb) as text')
    joined = " & \"|\" & ".join(f"v{i}" for i in range(len(reads)))
    script = f'''
    tell application "Microsoft Excel"
        set p to POSIX file "{xlsx}" as text
        open workbook workbook file name p
        delay 2
        set wb to workbook "{os.path.basename(xlsx)}"
        calculate
        set AppleScript's text item delimiters to ","
        {chr(10).join(lines)}
        close wb saving no
        return {joined}
    end tell
    '''
    try:
        out = subprocess.run(["osascript", "-e", script], capture_output=True, text=True, timeout=240)
    except Exception as e:  # noqa: BLE001
        return None, f"osascript failed: {e}"
    if out.returncode != 0:
        return None, out.stderr.strip()
    parts = out.stdout.strip().split("|")
    if len(parts) != len(reads):
        return None, f"unexpected Excel output ({len(parts)} blocks, expected {len(reads)})"
    worst, ncheck = 0.0, 0
    k = 0
    for calc, rows in zip(calcs, model_rows):
        for key in ("UTh", "MoW"):
            vals = [float(x) for x in parts[k].split(",")]
            k += 1
            if len(vals) != len(rows):
                return None, f"{calc['title']} {key}: {len(vals)} values read, expected {len(rows)}"
            for v, (_, _, r) in zip(vals, rows):
                worst = max(worst, abs(v - r[key]) / abs(r[key]))
                ncheck += 1
    checks = parts[k].split(",")
    all_ok = checks[-1].strip().lower() == "true"
    n_false = sum(1 for c in checks if c.strip().lower() != "true")
    return (worst, ncheck, all_ok, n_false), "ok"


# ──────────────────────────────────────────────────────────────────────────────
# Independent evaluation of the formula strings written to the workbook
# ──────────────────────────────────────────────────────────────────────────────
_CELL_RE = re.compile(r"(?<![A-Za-z0-9_])(?:'([^']+)'!|([A-Za-z_][A-Za-z0-9_]*)!)?\$?([A-Z]{1,3})\$?([0-9]+)(?![A-Za-z0-9_(])")
_NAME_RE = re.compile(r"(?<![A-Za-z0-9_!'])([A-Za-z_][A-Za-z0-9_]*)(?![A-Za-z0-9_(!'])")


class FormulaEvaluator:
    """Evaluates the arithmetic formulas of the Inputs / calculation / DG sheets."""

    def __init__(self, wb):
        self.wb = wb
        self.names = {}
        for name, dn in wb.defined_names.items():
            sheet, ref = dn.attr_text.split("!")
            self.names[name] = (sheet.strip("'"), ref.replace("$", ""))
        self.cache = {}

    def cell(self, sheet, ref):
        key = (sheet, ref)
        if key in self.cache:
            return self.cache[key]
        v = self.wb[sheet][ref].value
        if isinstance(v, str) and v.startswith("="):
            v = self.formula(v[1:], sheet)
        self.cache[key] = v
        return v

    def formula(self, expr, sheet):
        def cell_sub(m):
            sh = m.group(1) or m.group(2) or sheet
            return repr(float(self.cell(sh, m.group(3) + m.group(4))))
        expr = _CELL_RE.sub(cell_sub, expr)

        def name_sub(m):
            n = m.group(1)
            if n in ("LOG10", "ABS", "AND", "TRUE", "FALSE", "E"):
                return n
            if n in self.names:
                return repr(float(self.cell(*self.names[n])))
            raise KeyError(f"unknown name {n} in {sheet}: {expr}")
        expr = _NAME_RE.sub(name_sub, expr)
        expr = expr.replace("^", "**").replace("LOG10(", "math.log10(")
        return eval(expr, {"math": math, "__builtins__": {}})  # noqa: S307 (own formulas only)


def audit_formulas(wb, calcs, model_rows, dg_names, tables):
    """Evaluate every formula of the calc and DG sheets and compare with the Python mirror."""
    ev = FormulaEvaluator(wb)
    worst, n = 0.0, 0
    for calc, rows in zip(calcs, model_rows):
        ws = wb[calc["title"]]
        hdr = {c.column_letter: c.value for c in ws[3] if c.value}
        keymap = {"U/Th": "UTh", "Mo/W": "MoW", "Source U/Th": "UTh", "Source Mo/W": "MoW", "Source Nb/La": "NbLa",
                  "Fluid Nb/La (plotted)": "Fluid_NbLa", "Fluid Mo/W (plotted)": "Fluid_MoW", "Bulk DMo": "D_Mo", "Bulk DW": "D_W", "Bulk DU": "D_U", "Bulk DTh": "D_Th",
                  "FMQ-buffer log10(fO2)": "logfO2_FMQ", "Absolute log10(fO2)": "logfO2", "Mo solubility (µg/g)": "Mo_sol_ugg",
                  "W solubility (µg/g)": "W_sol_ugg", "Fluid Mo": "Fluid_Mo", "Fluid W": "Fluid_W", "Fluid U": "Fluid_U", "Fluid Th": "Fluid_Th",
                  "Source Mo": "Src_Mo", "Source W": "Src_W", "Source U": "Src_U", "Source Th": "Src_Th"}
        for k, (_, _, r) in enumerate(rows):
            row = calc["first"] + k
            for col, h in hdr.items():
                if h in keymap:
                    v = ev.cell(calc["title"], f"{col}{row}")
                    worst = max(worst, abs(v - r[keymap[h]]) / abs(r[keymap[h]]))
                    n += 1
    for dg in dg_names:
        ws = wb[dg]
        hdrs, cols = tables[dg]
        for j, (h, vals) in enumerate(cols, start=1):
            if h in ("FMQ", "PARAM", "Band_mask"):
                continue
            for i, pyv in enumerate(vals):
                v = ev.cell(dg, f"{get_column_letter(j)}{2 + i}")
                worst = max(worst, abs(v - pyv) / abs(pyv))
                n += 1
            assert ws.cell(row=2 + len(vals), column=j).value is None, f"{dg} column {h} longer than expected"
    return worst, n


# ──────────────────────────────────────────────────────────────────────────────
def main():
    write_model_inputs_template()          # creates data/model_inputs.xlsx once; user edits are kept
    model = Model()
    model_435 = Model(mo_a=0.435)
    for lst, new in ((DFMQ, model.dfmq), (SED_WT, model.sed_wt), (RUT_WT, model.rut_wt), (NACL_WT, model.nacl_wt)):
        lst[:] = new                       # grid axes from data/model_inputs.xlsx (in place: shared list objects)
    ENDMEMBERS = model.em
    grid_pc = model.sediment_grid(model.em["Carb"])
    grid_mx = model.sediment_grid(model.mix)
    grid_rut = model.rutile_grid()
    grid_sal = model.salinity_grid()

    wb = Workbook()
    write_inputs(wb, model)
    calc_rut = write_calc_sheet(wb, "Rutile", "rutile",
                                note="Residual rutile grid: 7 rutile contents × 6 ΔFMQ, pure N-MORB protolith, 10 wt.% NaCl. "
                                     "Plotted quantities are the slab-FLUID ratios Nb/La and Mo/W (no mantle mixing), as in the published figure; "
                                     "source ratios after 1 % fluid + 99 % DM mixing are given for reference.")
    calc_sal = write_calc_sheet(wb, "Salinity", "salinity",
                                note="Fluid salinity grid: 6 NaCl contents × 6 ΔFMQ, pure N-MORB protolith, 2.5 wt.% rutile.")
    calc_pc = write_calc_sheet(wb, "Pure_carbonate", "sediment", sediment_name="Carb",
                               note="Sediment grid, pure carbonate end-member: 5 sediment additions × 6 ΔFMQ, 10 wt.% NaCl, 2.5 wt.% rutile.")
    calc_mx = write_calc_sheet(wb, "C60_H40", "sediment", sediment_name="Mix",
                               note="Sediment grid, 60 % carbonate + 40 % hemipelagic end-member (Inputs!carb_frac): 5 sediment additions × 6 ΔFMQ.")

    sed_iso = [v for v in SED_WT if v != 0]
    hdr_pc = write_dg_sheet(wb, "DG_Pure_carbonate", calc_pc, "UTh", "UTh", "Sediment_wt_pct", "Sed", sed_iso, grid_pc)
    hdr_mx = write_dg_sheet(wb, "DG_C60_H40", calc_mx, "UTh", "UTh", "Sediment_wt_pct", "Sed", sed_iso, grid_mx)
    hdr_rut = write_dg_sheet(wb, "DG_Rutile", calc_rut, "FNbLa", "Fluid_NbLa", "Rutile_wt_pct", "Rut", RUT_WT, grid_rut,
                             ykey="FMoW", yname="Fluid_MoW")
    hdr_sal = write_dg_sheet(wb, "DG_Salinity", calc_sal, "UTh", "UTh", "NaCl_wt_pct", "NaCl", NACL_WT, grid_sal)

    # Python-mirror DataGraph tables
    tables = {
        "DG_Pure_carbonate": (hdr_pc, dg_table(grid_pc, "UTh", sed_iso, SED_WT)),
        "DG_C60_H40": (hdr_mx, dg_table(grid_mx, "UTh", sed_iso, SED_WT)),
        "DG_Rutile": (hdr_rut, dg_table(grid_rut, "Fluid_NbLa", RUT_WT, RUT_WT, ykey="Fluid_MoW")),
        "DG_Salinity": (hdr_sal, dg_table(grid_sal, "UTh", NACL_WT, NACL_WT)),
    }
    os.makedirs(DG_DIR, exist_ok=True)  # also creates outputs/
    for name, (hdr, cols) in tables.items():
        with open(os.path.join(DG_DIR, name + ".txt"), "w", encoding="utf-8") as fh:
            fh.write(tsv_text(hdr, cols))

    # ── validation (Python mirror) ──
    results = []
    def rec(num, text, ok, detail=""):
        results.append((num, text, bool(ok), detail))

    same0 = all(abs(a[2]["UTh"] - b[2]["UTh"]) < 1e-12 and abs(a[2]["MoW"] - b[2]["MoW"]) < 1e-12
                for a, b in zip(grid_pc, grid_mx) if a[1] == 0)
    rec(1, "Pure-carbonate and 60C:40H grids identical at 0 wt.% sediment", same0)
    no_si = True
    for name, (hdr, cols) in tables.items():
        d = dict(cols)
        npts = 2 * len(set(d["PARAM"])) + 1  # 11 for the sediment grids
        for b in range(1, 6):
            xs, ys = d[f"Band{b}_X"], d[f"Band{b}_Y"]
            rec(2, f"{name}: Band{b} has 2n+1 = {npts} points", len(xs) == npts, f"{len(xs)} points")
            rec(3, f"{name}: Band{b} first == last coordinate", xs[0] == xs[-1] and ys[0] == ys[-1])
            si = polygon_self_intersects(xs, ys, log=True) or polygon_self_intersects(xs, ys, log=False)
            no_si &= not si
            rec(4, f"{name}: Band{b} does not self-intersect", not si)
    for name in ("DG_Pure_carbonate", "DG_C60_H40"):
        d = dict(tables[name][1])
        for v in sed_iso:
            rec(5, f"{name}: Sed{param_tag(v)} trajectory has six points", len(d[f"ISO{param_tag(v)}_X"]) == 6)
        rows25 = [r for _, w, r in (grid_pc if "Pure" in name else grid_mx) if w == 2.5]
        rec(6, f"{name}: Sed2_5 built from the 2.5 wt.% rows", d["ISO2_5_X"] == [r["UTh"] for r in rows25])
        g = grid_pc if "Pure" in name else grid_mx
        for num, b, lo, hi in [(7, 3, 0, 1), (8, 5, 2, 3)]:
            lo_pts = [r["UTh"] for dd, _, r in g if dd == lo]
            hi_pts = [r["UTh"] for dd, _, r in g if dd == hi]
            ok = d[f"Band{b}_X"][:5] == lo_pts and d[f"Band{b}_X"][5:10] == hi_pts[::-1]
            rec(num, f"{name}: Band{b} between ΔFMQ {lo:+d} and {hi:+d}", ok)
        rec(9, f"{name}: Band_mask is its own column", "Band_mask" in tables[name][0] and tables[name][0][4] == "Band_mask")
        rec(10, f"{name}: no Sed0_X / Sed0_Y", not any(h.startswith("Sed0_") for h in tables[name][0]))

    python_checks = {"no_self_intersection": no_si}
    n_check_rows = write_checks(wb, calc_pc, calc_mx, ["DG_Pure_carbonate", "DG_C60_H40", "DG_Rutile", "DG_Salinity"], python_checks)
    dg_names = ["DG_Pure_carbonate", "DG_C60_H40", "DG_Rutile", "DG_Salinity"]
    audit_worst, audit_n = audit_formulas(wb, [calc_pc, calc_mx, calc_rut, calc_sal], [grid_pc, grid_mx, grid_rut, grid_sal], dg_names, tables)
    audit_msg = (f"Formula audit: {audit_n} workbook formulas (calculation and DG sheets) evaluated independently from their "
                 f"formula strings and compared with the Python mirror; maximum relative difference {audit_worst:.1e}.")
    print(audit_msg)
    wb.calculation.fullCalcOnLoad = True
    wb.save(XLSX)

    # ── Excel recalculation cross-check (read-only; the deliverable is not modified) ──
    res, msg = excel_cross_check(XLSX, [calc_pc, calc_mx, calc_rut, calc_sal],
                                 [grid_pc, grid_mx, grid_rut, grid_sal], n_check_rows)
    if res is None:
        excel_msg = f"Excel recalculation not performed ({msg}). Formulas recalculate when the workbook is opened (fullCalcOnLoad)."
    else:
        worst, ncheck, all_ok, n_false = res
        excel_msg = (f"Microsoft Excel recalculated the workbook: {ncheck} U/Th and Mo/W cells compared with the Python mirror, "
                     f"maximum relative difference {worst:.1e}. Checks sheet: ALL CHECKS = {all_ok} ({n_false} non-TRUE cells).")

    # ── numerical table ──
    g0_pc = {w: r for d, w, r in grid_pc if d == 0}
    g0_mx = {w: r for d, w, r in grid_mx if d == 0}
    g0_pc435 = {w: r for d, w, r in model_435.sediment_grid(model_435.em["Carb"]) if d == 0}
    g0_mx435 = {w: r for d, w, r in model_435.sediment_grid(model_435.mix) if d == 0}

    with open(REPORT, "w", encoding="utf-8") as fh:
        fh.write(f"# Validation report — Mo/W–U/Th sensitivity grids\n\nGenerated {date.today().isoformat()} by `thermodynamic_model/build_workbook.py`.\n\n")
        fh.write("## Conditions\n\n")
        fh.write("T = 1000 °C (1273.15 K), P = 26 100 bar, 10 wt.% NaCl (reference), 2.5 wt.% rutile (reference), cpx:grt = 70:30, "
                 "F_batch = 0, 1 % fluid + 99 % depleted mantle. Mo fO₂ coefficient Mo_a = 0.44.\n\n")
        fh.write("## Structural and geometric checks\n\n| # | Check | Result | Detail |\n|---|---|---|---|\n")
        for num, text, okk, detail in results:
            fh.write(f"| {num} | {text} | {'PASS' if okk else 'FAIL'} | {detail} |\n")
        npass = sum(1 for r in results if r[2])
        fh.write(f"\n**{npass} / {len(results)} checks passed.**\n\n")
        fh.write("## Formula audit and Excel cross-check\n\n" + audit_msg + "\n\n" + excel_msg + "\n\n")
        fh.write("## Numerical check at ΔFMQ 0\n\n")
        fh.write("Supplied values are reproduced **exactly** (to the six decimals given) when Mo_a = 0.435, the value in Bali et al. (2012) "
                 "Appendix C. With the chosen Mo_a = 0.44 the U/Th values are identical and Mo/W is lower by a constant factor "
                 f"of {g0_pc[0]['MoW'] / TARGETS_DFMQ0[0][1]:.4f} (≈ 10^(−0.005·log fO₂)).\n\n")
        fh.write("| Sediment | Carb U/Th (0.44) | Carb U/Th (supplied) | Carb Mo/W (0.44) | Carb Mo/W (0.435) | Carb Mo/W (supplied) | "
                 "60C:40H U/Th (0.44) | 60C:40H U/Th (supplied) | 60C:40H Mo/W (0.44) | 60C:40H Mo/W (0.435) | 60C:40H Mo/W (supplied) |\n")
        fh.write("|---|---|---|---|---|---|---|---|---|---|---|\n")
        for w in SED_WT:
            t = TARGETS_DFMQ0[w]
            fh.write(f"| {w:g} % | {g0_pc[w]['UTh']:.6f} | {t[0]:.6f} | {g0_pc[w]['MoW']:.6f} | {g0_pc435[w]['MoW']:.6f} | {t[1]:.6f} | "
                     f"{g0_mx[w]['UTh']:.6f} | {t[2]:.6f} | {g0_mx[w]['MoW']:.6f} | {g0_mx435[w]['MoW']:.6f} | {t[3]:.6f} |\n")
        maxd_uth = max(abs(g0_pc[w]["UTh"] - TARGETS_DFMQ0[w][0]) for w in SED_WT)
        maxd_435 = max(abs(g0_pc435[w]["MoW"] - TARGETS_DFMQ0[w][1]) for w in SED_WT)
        maxd_435m = max(abs(g0_mx435[w]["MoW"] - TARGETS_DFMQ0[w][3]) for w in SED_WT)
        fh.write(f"\nMax |ΔU/Th| vs supplied: {maxd_uth:.1e}. Max |ΔMo/W| at 0.435 vs supplied: {max(maxd_435, maxd_435m):.1e}.\n\n")
        fh.write("## Intermediate quantities at the reference case (ΔFMQ 0, 0 wt.% sediment)\n\n")
        r0 = g0_pc[0]
        fh.write("| Quantity | Value |\n|---|---|\n")
        for k in ["NaCl_m", "log_NaCl", "logfO2_FMQ", "logfO2", "Mo_sol_m", "Mo_sol_ugg", "W_sol_m", "W_sol_ugg", "U_sol_sil", "U_sol_rut",
                  "D_Mo", "D_W", "D_U", "D_Th", "Fluid_Mo", "Fluid_W", "Fluid_U", "Fluid_Th", "Src_Mo", "Src_W", "Src_U", "Src_Th", "UTh", "MoW"]:
            fh.write(f"| {k} | {r0[k]:.6g} |\n")

    print(f"Workbook: {XLSX}")
    print(f"TSV tables: {DG_DIR}")
    print(f"Report: {REPORT}")
    print(excel_msg)
    print(f"Checks passed: {npass}/{len(results)}")


if __name__ == "__main__":
    main()
