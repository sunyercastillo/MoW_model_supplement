# Validation report — Mo/W–U/Th sensitivity grids

Generated 2026-10-30 by `thermodynamic_model/build_workbook.py`.

## Conditions

T = 1000 °C (1273.15 K), P = 26 100 bar, 10 wt.% NaCl (reference), 2.5 wt.% rutile (reference), cpx:grt = 70:30, F_batch = 0, 1 % fluid + 99 % depleted mantle. Mo fO₂ coefficient Mo_a = 0.44.

## Structural and geometric checks

| # | Check | Result | Detail |
|---|---|---|---|
| 1 | Pure-carbonate and 60C:40H grids identical at 0 wt.% sediment | PASS |  |
| 2 | DG_Pure_carbonate: Band1 has 2n+1 = 11 points | PASS | 11 points |
| 3 | DG_Pure_carbonate: Band1 first == last coordinate | PASS |  |
| 4 | DG_Pure_carbonate: Band1 does not self-intersect | PASS |  |
| 2 | DG_Pure_carbonate: Band2 has 2n+1 = 11 points | PASS | 11 points |
| 3 | DG_Pure_carbonate: Band2 first == last coordinate | PASS |  |
| 4 | DG_Pure_carbonate: Band2 does not self-intersect | PASS |  |
| 2 | DG_Pure_carbonate: Band3 has 2n+1 = 11 points | PASS | 11 points |
| 3 | DG_Pure_carbonate: Band3 first == last coordinate | PASS |  |
| 4 | DG_Pure_carbonate: Band3 does not self-intersect | PASS |  |
| 2 | DG_Pure_carbonate: Band4 has 2n+1 = 11 points | PASS | 11 points |
| 3 | DG_Pure_carbonate: Band4 first == last coordinate | PASS |  |
| 4 | DG_Pure_carbonate: Band4 does not self-intersect | PASS |  |
| 2 | DG_Pure_carbonate: Band5 has 2n+1 = 11 points | PASS | 11 points |
| 3 | DG_Pure_carbonate: Band5 first == last coordinate | PASS |  |
| 4 | DG_Pure_carbonate: Band5 does not self-intersect | PASS |  |
| 2 | DG_C60_H40: Band1 has 2n+1 = 11 points | PASS | 11 points |
| 3 | DG_C60_H40: Band1 first == last coordinate | PASS |  |
| 4 | DG_C60_H40: Band1 does not self-intersect | PASS |  |
| 2 | DG_C60_H40: Band2 has 2n+1 = 11 points | PASS | 11 points |
| 3 | DG_C60_H40: Band2 first == last coordinate | PASS |  |
| 4 | DG_C60_H40: Band2 does not self-intersect | PASS |  |
| 2 | DG_C60_H40: Band3 has 2n+1 = 11 points | PASS | 11 points |
| 3 | DG_C60_H40: Band3 first == last coordinate | PASS |  |
| 4 | DG_C60_H40: Band3 does not self-intersect | PASS |  |
| 2 | DG_C60_H40: Band4 has 2n+1 = 11 points | PASS | 11 points |
| 3 | DG_C60_H40: Band4 first == last coordinate | PASS |  |
| 4 | DG_C60_H40: Band4 does not self-intersect | PASS |  |
| 2 | DG_C60_H40: Band5 has 2n+1 = 11 points | PASS | 11 points |
| 3 | DG_C60_H40: Band5 first == last coordinate | PASS |  |
| 4 | DG_C60_H40: Band5 does not self-intersect | PASS |  |
| 2 | DG_Rutile: Band1 has 2n+1 = 15 points | PASS | 15 points |
| 3 | DG_Rutile: Band1 first == last coordinate | PASS |  |
| 4 | DG_Rutile: Band1 does not self-intersect | PASS |  |
| 2 | DG_Rutile: Band2 has 2n+1 = 15 points | PASS | 15 points |
| 3 | DG_Rutile: Band2 first == last coordinate | PASS |  |
| 4 | DG_Rutile: Band2 does not self-intersect | PASS |  |
| 2 | DG_Rutile: Band3 has 2n+1 = 15 points | PASS | 15 points |
| 3 | DG_Rutile: Band3 first == last coordinate | PASS |  |
| 4 | DG_Rutile: Band3 does not self-intersect | PASS |  |
| 2 | DG_Rutile: Band4 has 2n+1 = 15 points | PASS | 15 points |
| 3 | DG_Rutile: Band4 first == last coordinate | PASS |  |
| 4 | DG_Rutile: Band4 does not self-intersect | PASS |  |
| 2 | DG_Rutile: Band5 has 2n+1 = 15 points | PASS | 15 points |
| 3 | DG_Rutile: Band5 first == last coordinate | PASS |  |
| 4 | DG_Rutile: Band5 does not self-intersect | PASS |  |
| 2 | DG_Salinity: Band1 has 2n+1 = 13 points | PASS | 13 points |
| 3 | DG_Salinity: Band1 first == last coordinate | PASS |  |
| 4 | DG_Salinity: Band1 does not self-intersect | PASS |  |
| 2 | DG_Salinity: Band2 has 2n+1 = 13 points | PASS | 13 points |
| 3 | DG_Salinity: Band2 first == last coordinate | PASS |  |
| 4 | DG_Salinity: Band2 does not self-intersect | PASS |  |
| 2 | DG_Salinity: Band3 has 2n+1 = 13 points | PASS | 13 points |
| 3 | DG_Salinity: Band3 first == last coordinate | PASS |  |
| 4 | DG_Salinity: Band3 does not self-intersect | PASS |  |
| 2 | DG_Salinity: Band4 has 2n+1 = 13 points | PASS | 13 points |
| 3 | DG_Salinity: Band4 first == last coordinate | PASS |  |
| 4 | DG_Salinity: Band4 does not self-intersect | PASS |  |
| 2 | DG_Salinity: Band5 has 2n+1 = 13 points | PASS | 13 points |
| 3 | DG_Salinity: Band5 first == last coordinate | PASS |  |
| 4 | DG_Salinity: Band5 does not self-intersect | PASS |  |
| 5 | DG_Pure_carbonate: Sed1 trajectory has six points | PASS |  |
| 5 | DG_Pure_carbonate: Sed2_5 trajectory has six points | PASS |  |
| 5 | DG_Pure_carbonate: Sed5 trajectory has six points | PASS |  |
| 5 | DG_Pure_carbonate: Sed10 trajectory has six points | PASS |  |
| 6 | DG_Pure_carbonate: Sed2_5 built from the 2.5 wt.% rows | PASS |  |
| 7 | DG_Pure_carbonate: Band3 between ΔFMQ +0 and +1 | PASS |  |
| 8 | DG_Pure_carbonate: Band5 between ΔFMQ +2 and +3 | PASS |  |
| 9 | DG_Pure_carbonate: Band_mask is its own column | PASS |  |
| 10 | DG_Pure_carbonate: no Sed0_X / Sed0_Y | PASS |  |
| 5 | DG_C60_H40: Sed1 trajectory has six points | PASS |  |
| 5 | DG_C60_H40: Sed2_5 trajectory has six points | PASS |  |
| 5 | DG_C60_H40: Sed5 trajectory has six points | PASS |  |
| 5 | DG_C60_H40: Sed10 trajectory has six points | PASS |  |
| 6 | DG_C60_H40: Sed2_5 built from the 2.5 wt.% rows | PASS |  |
| 7 | DG_C60_H40: Band3 between ΔFMQ +0 and +1 | PASS |  |
| 8 | DG_C60_H40: Band5 between ΔFMQ +2 and +3 | PASS |  |
| 9 | DG_C60_H40: Band_mask is its own column | PASS |  |
| 10 | DG_C60_H40: no Sed0_X / Sed0_Y | PASS |  |

**79 / 79 checks passed.**

## Formula audit and Excel cross-check

Formula audit: 3638 workbook formulas (calculation and DG sheets) evaluated independently from their formula strings and compared with the Python mirror; maximum relative difference 4.9e-16.

Excel recalculation not performed (261:306: execution error: Microsoft Excel got an error: The object you are trying to access does not exist (-1728)). Formulas recalculate when the workbook is opened (fullCalcOnLoad).

## Numerical check at ΔFMQ 0

Supplied values are reproduced **exactly** (to the six decimals given) when Mo_a = 0.435, the value in Bali et al. (2012) Appendix C. With the chosen Mo_a = 0.44 the U/Th values are identical and Mo/W is lower by a constant factor of 0.9250 (≈ 10^(−0.005·log fO₂)).

| Sediment | Carb U/Th (0.44) | Carb U/Th (supplied) | Carb Mo/W (0.44) | Carb Mo/W (0.435) | Carb Mo/W (supplied) | 60C:40H U/Th (0.44) | 60C:40H U/Th (supplied) | 60C:40H Mo/W (0.44) | 60C:40H Mo/W (0.435) | 60C:40H Mo/W (supplied) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 % | 0.754334 | 0.754334 | 7.558422 | 8.170971 | 8.170971 | 0.754334 | 0.754334 | 7.558422 | 8.170971 | 8.170971 |
| 1 % | 0.792138 | 0.792138 | 7.405431 | 8.007876 | 8.007876 | 0.871773 | 0.871773 | 7.221007 | 7.814167 | 7.814167 |
| 2.5 % | 0.848811 | 0.848811 | 7.193770 | 7.782239 | 7.782239 | 1.045135 | 1.045135 | 6.812702 | 7.382399 | 7.382399 |
| 5 % | 0.943183 | 0.943183 | 6.881691 | 7.449550 | 7.449550 | 1.326850 | 1.326850 | 6.308544 | 6.849271 | 6.849271 |
| 10 % | 1.131609 | 1.131609 | 6.375501 | 6.909933 | 6.909933 | 1.864670 | 1.864670 | 5.670292 | 6.174343 | 6.174343 |

Max |ΔU/Th| vs supplied: 4.1e-07. Max |ΔMo/W| at 0.435 vs supplied: 4.2e-07.

## Intermediate quantities at the reference case (ΔFMQ 0, 0 wt.% sediment)

| Quantity | Value |
|---|---|
| NaCl_m | 1.90129 |
| log_NaCl | 0.279047 |
| logfO2_FMQ | -8.72202 |
| logfO2 | -8.72202 |
| Mo_sol_m | 0.46312 |
| Mo_sol_ugg | 44436.3 |
| W_sol_m | 1.27759 |
| W_sol_ugg | 234872 |
| U_sol_sil | 363.494 |
| U_sol_rut | 47.304 |
| D_Mo | 0.0500167 |
| D_W | 0.0314393 |
| D_U | 0.10252 |
| D_Th | 0.9931 |
| Fluid_Mo | 8.15727 |
| Fluid_W | 1.20868 |
| Fluid_U | 0.692548 |
| Fluid_Th | 0.188299 |
| Src_Mo | 0.106323 |
| Src_W | 0.0140668 |
| Src_U | 0.0118755 |
| Src_Th | 0.015743 |
| UTh | 0.754334 |
| MoW | 7.55842 |
