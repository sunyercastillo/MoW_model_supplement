# Audit of inputs and equations — Mo/W–U/Th sensitivity grids

Scope: the forward model of Bali et al. (2011, 2012) as implemented in the
project workbook `Bali model.xlsx` (sheet `Calculations`), in the Python engine
`CORE_scientific_calculations_engine.py`, and as published in Bali et al. (2012)
Appendix C (`1-s2.0-S0012821X12004074-mmc3.doc`). Date: 2026-10-30.

## 1. Equations used (identical in workbook and engine)

| Step | Expression | Source |
|---|---|---|
| NaCl molality | m = 1000·w / (58.44·(100 − w)), w in wt.% | workbook `Calculations!C17:C24` |
| FMQ buffer | log10 fO₂(FMQ) = −25096.3/T + 8.735 + 0.11·(P − 1)/T, T in K, **P in bar** | workbook `D7`; engine `compute_base_fmq_log10` |
| Absolute fO₂ | log10 fO₂ = log10 fO₂(FMQ) + ΔFMQ | workbook `E7` |
| Mo solubility | log10 m_Mo = a·log fO₂ + 0.42·log10 m_NaCl − 1.8·(1000/T) + 4.8; c_Mo (µg/g) = m_Mo·95.95·1000 | workbook `K3:M3`; App. C |
| W solubility | log10 m_W = 0.07·log fO₂ − 4.7236·(1000/T) + 4.4271; c_W (µg/g) = m_W·183.84·1000 | workbook `K5:L5`; App. C |
| U solubility | log10 c_U = 2.681 + 0.1433·log fO₂ + 0.594·m_NaCl (silicates); intercept 1.7954 for rutile | workbook `K6:M6`; App. C; Bali et al. (2011) |
| Mineral/fluid D (Mo, W, U) | D = c_sat / c_solubility with c_sat: Mo 40 (cpx), 12 (grt), 87 670 (rut); W 60 (cpx), 12 (grt); U 11 (cpx), 40 (grt), 94 (rut) | App. C |
| Constant D | W rut 1.25; Th 1.19 / 0.61 / 0.10; Nb 0.172 / 0.204 / 200; La 1.429 / 0.204 / 1.25 (cpx / grt / rut) | App. C |
| Bulk D | D_bulk = x_cpx·D_cpx + x_grt·D_grt + x_rut·D_rut | workbook `N3:N9` |
| Modes | x_rut = wt.%/100; x_cpx = 0.70·(1 − x_rut); x_grt = 0.30·(1 − x_rut) | workbook `FMQ!G64:N67` |
| Fluid | C_fluid = C₀ / (F + (1 − F)·D_bulk), F = 0 | workbook `G13:M13`, `N13 = 0` |
| Source | C_source = f·C_fluid + (1 − f)·C_DM, f = 0.01 | workbook `F24:M55` |
| Protolith with sediment | C₀ = (1 − s)·C_NMORB + s·C_sed | `sediment_endmembers.py`; Supplementary Methods Eq. 9.1 |
| Mixed sediment | C_60C40H = 0.6·C_carb + 0.4·C_hemi | task specification |

## 2. Inconsistencies found

### 2.1 Mo fO₂ coefficient: 0.435 (Appendix C) vs 0.44 (workbook D cells, engine, methods text)

Bali et al. (2012) Appendix C gives all three Mo partition coefficients as
`c_sat / 10^(0.435·log fO₂ + 0.42·log NaCl − 1.8·1000/T + 4.8)`. The project
workbook uses **0.44** in the partition-coefficient cells (`Calculations!K3:M3`)
but **0.435** in auxiliary cells (`Calculations!T3`, `FMQ!W2`, `Sheet3!D2`);
the Python engine, the Supplementary Methods (Eq. 3) and all project notes use
0.44.

The check values supplied with the task (e.g. U/Th 0.754334, Mo/W 8.170971 at
ΔFMQ 0, 0 wt.% sediment) are reproduced to six decimals **only with 0.435**.
With 0.44, U/Th is unchanged and Mo/W is lower by the constant factor
10^(−0.005·log fO₂) = 0.925 at ΔFMQ 0 (7.558 instead of 8.171).

**Decision (author):** keep 0.44 in the workbook and the engine. The
coefficient is a single named cell (`Inputs!Mo_a`); setting it to 0.435
reproduces the supplied table exactly. The validation report lists both.

### 2.2 End-member precision

| Element | Table B1 (used) | Workbook / Appendix C |
|---|---|---|
| N-MORB U | 0.071 | 0.0711 |
| N-MORB Th | 0.187 | 0.1871 |
| DM W | 0.002 | 0.0024 |
| DM U | 0.005 | 0.0047 |
| DM Th | 0.014 | 0.0137 |

Effect at ΔFMQ 0, 0 wt.% sediment, Mo_a = 0.44: Table B1 values give U/Th
0.7543 and Mo/W 7.558; full-precision values give U/Th 0.7502 and Mo/W 7.352
(−2.7 %). The DM tungsten value (0.002 vs 0.0024) is responsible for nearly all
of the Mo/W difference, because the 99 % mantle term dominates the source W
budget. **Decision (author):** Table B1 values.

### 2.3 Sediment end-members in the Supplementary Methods text

`Supplementary_Methods.md` §9 still lists older sediment values (carbonate Mo
1.09, W 0.19, U 0.57; hemipelagic Mo 3.96, W 1.22, U 3.72; no Th). The code
(`sediment_endmembers.py`) and Table B1 use the EA-1 unit means: carbonate
1.087 / 0.187 / 0.684 / 0.222 and hemipelagic 5.626 / 0.874 / 4.117 / 2.656
(Mo / W / U / Th). The grids use Table B1. The methods text must be updated.

### 2.4 Legacy N-MORB dictionary in the engine

`BaliConstants.morb_concentrations` (Mo 0.19, W 0.040, U 0.047, Th 0.087) is a
legacy default that is overridden whenever the input workbook supplies C₀. It is
not used by the grids. It should be removed or replaced by the Table B1 values
before publication to avoid confusion.

### 2.5 Panel A (residual rutile) uses the fluid composition

The published figure plots panel A with the slab-fluid ratios (Nb/La = 5.26 at
0 wt.% rutile, 0.185 at 2.5 wt.%, 0.095 at 5 wt.%), whereas panels B–D plot the
metasomatised source (1 % fluid + 99 % DM). With source mixing the rutile grid
would be compressed to Nb/La 0.79–1.49 because DM Nb/La (0.90) dominates.
**Decision (author):** keep panel A as the fluid composition and label it as such
(`Fluid_NbLa`, `Fluid_MoW`); the `Rutile` sheet also lists the source ratios.

### 2.6 Minor

- Appendix C lists subducted-sediment La as 28.8 µg/g; the workbook has 29.8.
  Not used by the grids.
- The workbook's `Calculations!T3` and the auxiliary Mo cell in `FMQ!W2` use
  0.435 while the D cells use 0.44 (see 2.1).

## 3. Values used in the grids

Conditions: T = 1000 °C (1273.15 K); P = 26 100 bar; NaCl = 10 wt.%
(m = 1.9013 mol/kg); rutile 2.5 wt.% (x_cpx = 0.6825, x_grt = 0.2925,
x_rut = 0.025); ΔFMQ = −2, −1, 0, +1, +2, +3; F = 0; f = 0.01.

End-members (µg/g): N-MORB Mo 0.408, W 0.038, U 0.071, Th 0.187, Nb 3.507,
La 3.895; DM Mo 0.025, W 0.002, U 0.005, Th 0.014, Nb 0.21, La 0.234;
carbonate Mo 1.087, W 0.187, U 0.684, Th 0.222; hemipelagic Mo 5.626, W 0.874,
U 4.117, Th 2.656; 60C:40H Mo 2.9026, W 0.4618, U 2.0572, Th 1.1956.

Reference-case intermediates (ΔFMQ 0, 0 wt.% sediment): log fO₂ = −8.7220;
c_Mo = 44 436 µg/g; c_W = 234 872 µg/g; c_U(sil) = 363.5, c_U(rut) = 47.30
µg/g; D_Mo = 0.05002, D_W = 0.03144, D_U = 0.10252, D_Th = 0.9931; fluid
Mo 8.157, W 1.2087, U 0.6925, Th 0.1883; source U/Th 0.7543, Mo/W 7.558.

## 4. Differences from the earlier Bali workbook and grid

1. The sediment protolith, the 60C:40H mixture and the four grids are new; the
   earlier workbook computed single cases and a fluid–DM mixing table.
2. All constants are named cells in `Inputs`; the earlier workbook mixed
   constants into formulas (and used two different Mo coefficients).
3. The batch-fluid expression is kept with F as an explicit input (F = 0).
4. Panel A keeps the fluid composition of the published figure; this is now
   stated explicitly rather than implied.
5. Check values supplied for ΔFMQ 0 are reproduced with Mo_a = 0.435 and
   differ by a constant factor 0.925 in Mo/W with the retained 0.44.
