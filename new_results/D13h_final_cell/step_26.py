result = r"""# D13h: the final analytic cell N >= 10^10 and the [100,400] tightening — THEOREM for all N >= 100

**Status: COMPLETE.** All three ranges are certified by arb enclosures (python-flint, ctx.dps=70).

## 1. THEOREM (PROVED (arb), modulo §4.3's classical Binet formula)

For every integer N >= 100,

  |M_k^main + (sqrt(pi)/N) e^{-sqrt(pi N)} cos sqrt(pi N)|  <=  C * N^{-3/2} e^{-sqrt(pi N)}

with

  **C = 617117.6**

= max of the three certified suprema below (each an arb upper bound):

| range | cover | sup C | where |
|---|---|---|---|
| [100, 400] | 29 ratio-1.05 N-cells, `step_38.certified_y38` (J=24, nsub=4) | **194.3809** | cell [373.3456, 392.0129] |
| [400, 10^10] | 43 exactly-contiguous cells (planner's adaptive sweep, ratio 1.5 + split) | **617117.5536** | cell [6.6344e9, 9.9516e9] |
| [10^10, oo) | ONE analytic cell y in [0, y_J], y_J = sqrt(pi/2)*1e-5 (this work, `step_10_J192.py`) | **49.5626** | — |

The sweep constant grows like a power of N on wide boxes (tail factor e^{-c2R d^2} = y^{2r}, r = c2R/c2 < 1
over a ratio-1.5 box); it is a certified but lossy cover. The final cell avoids the loss by exact
identities, hence the much smaller constant on the infinite range.

## 2. The final analytic cell y in [0, y_J] (Task 1) — PROVED (arb)

`step_10.py` (J=48, C=53.6157) and `step_10_J192.py` (J=192, **C=49.5626**); J only refines, every J
certifies. Every piece of step_38's `cell()` is written as y^3 * (regular factor) and every regular
factor is evaluated on the box y in [0, y_J], sigma in [0,1], rho = 1+y^2 sigma in [1, 1+y_J^2],
t = rho y in [0, T], T = rho_hat y_J = 1.25334e-5, with **no division by y** (powers of y are tracked
symbolically; w^2 = log N + c_w = a_ + 2L, a_ = log(pi/2)+c_w, L = log(1/y), w_J = 4.901617).

### 2.1 Factored forms used (all exact algebra)
- Phi3b = y^2 Ph3reg, Ph3reg = pi rho^4 (rho^2 y^2 + 3 eta rho y + 3 eta^2)/(eta^3 (eta^2+rt2 eta rho y+rho^2 y^2)^{3/2})
- Phi4b = y^3 Ph4reg (same denominator structure)
- G1 = y G1reg, G1reg = 6 rho^2 y/eta^2 + (1.5 rt2 + 1) rho/eta;  G2 = y^2 G2reg, G2reg = 12 rt2 rho^3 y/eta^3 + 4 rho^2/eta^2
- c2 = y c2reg, c2reg = (pi rho^3/4) numR/(1+rt2 t+t^2)^2, numR = t^3+2 rt2 t^2+4 t+rt2; c2reg box = [1.1107 +/- 6.0e-5]
- om = y omreg, omreg = (pi/(2 rt2)) Omn/den1, box = [1.5708 +/- 5.9e-5]; mS = y mSreg, mSreg_up = pi rho^3 rt2/12
- **g0 via h1-Binet** (§4.3): |g0| <= t^{5/2} (sqrt pi/4) 2^{3/2} e^{ReE},
  ReE <= (t/24)(-rt2 + 459 t^2/(1-3T)) <= 0 on the box (asserted), so P0t <= t^{5/2} P0reg, P0reg = (sqrt pi/4) 2^{3/2} e^{y_J mSreg_up}.
- v-moments: v = u/x0 = sqrt(y) z, K2 = c2reg/(rho^2 y) = K2reg/y, K2reg = c2reg/rho^2; the Gaussian is
  e^{-K2reg z^2} exactly; z <= zmax = w rho/sqrt(2 c2reg); zmax_lo = 3.2886 > zstar = 3 (asserted).
- tails: d = w/sqrt(2 c2) (< x0/3, asserted), c2 d^2 = w^2/2 exactly, e^{-c2 d^2} = y sqrt(2/pi) e^{-c_w/2}
  exactly; e^{(c2-c2R)d^2} <= E1R = 1.87582 using |c2(1)-c2(eta)| <= y Lip_eta sqrt(y) zmax,
  Lip_eta = 4.9953 (interval quotient bound for |d/deta curv2reg|, true value ~3.33), and
  sqrt(y) w^3 <= sqrt(y_J) w_J^3 (monotonicity M2).
- lamR/lamL >= sqrt(y) w_J lamreg, lamreg = (pi rho^2/2) g1/sqrt(2 c2reg) - a y_J^{3/2} |Q1|/w_J:
  g1R_lo = 1.31901, g1L_lo = 1.36506 from the explicit -Gh' = (rho/rt2)(2 eta^3+4 rt2 eta^2 t+4 eta t^2+ rt2 t^3)/(eta(eta^2+rt2 eta t+t^2))^2;
  lamRreg = 1.39009 > 0, lamLreg = 1.43862 > 0 (asserted).
- TLf: ReX2 = th (M1-L2), MmL_up = -0.249982 < 0 (asserted) => TLf*N^{3/2} <= (pi/2)^{3/2} G1h e^{-K_J/y_J} y_J^{-3}/lam2
  with K_J = pi|Mml_up|/rt2 = 0.5553, y_J < K_J/3 (asserted, monotonicity M3): value 2.8e-19229.
- TL0: decreasing for N >= N_J (LpN < 0 asserted, L''<0 for N>=6 analytically): value 9.9e-2666376435.
- D: z-factorisation z = e^{i pi/4} rho zt3 e^E e^{mS} e^{om^2 invp2/2} e^{i phi} (verified EXACT, ratio 1 +/- 1.3e-57),
  |z-1| <= exp(sum d_k y) - 1 <= y (sum d_k) e^{y_J sum d_k}, d1..d6 = 1.3e-5, 1.22142, 1.50002, 0.37024, 0.39272, 1.66609;
  D <= (pi/rt2)(sum d_k) e^{y_J sum d_k} = Dreg = 11.4423.

### 2.2 Pieces (C-parts = (pi/2)^{3/2} * regular factor, i.e. piece*N^{3/2} bounds; C = (2/pi) sum + Dreg)

| piece | J=48 | J=192 | validation |
|---|---|---|---|
| W1 | 15.3849 | 15.3849 | EXACT match to step_38 W1*N^{1.5} at N=1e9 (15.383) |
| W2 | 9.3001 | 8.1733 | true integral 7.70 <= ours <= step_38's 11.15 |
| W3 | 38.3522 | 33.1124 | true integral 28.85 <= ours <= step_38's 46.81 |
| TM | 0.2026 | 0.2026 | step_38: 0.255 |
| TR | 0.4649 | 0.4649 | step_38: 0.531 |
| TLm | 2.5412 | 2.5412 | step_38: 1.604 (ours uses crude 2^{19/4} GLm2reg) |
| TLf | 2.8e-19229 | 2.8e-19229 | super-exponential |
| TL0 | 9.9e-2666376435 | 9.9e-2666376435 | super-exponential |
| D | 11.4423 | 11.4423 | step_38: ~2.44 (ours uses crude |E| bound) |
| **C** | **53.6157** | **49.5626** | pointwise limit C(1e9) = 51.16 decreasing |

### 2.3 Asserts (all pass, printed by step_10.py)
Omn > 0; 459 T^2/(1-3T) < rt2 (=> ReE <= 0); y_J(a_+2L_J) < 2 c2reg_lo/(9 rho_hat^2) (d < x0/3);
zstar <= zmax_lo; lamRreg, lamLreg > 0; denreg_lo = 10.29 > 0 (TM); lam2 = 3.332 > 0;
MmL_up < 0; y_J < K_J/3; LpN < 0 (TL0); a_+2L_J = 23.92 > 6 (monotonicity windows M1, M2).

### 2.4 The one classical input (label: PROVED modulo Binet; NUMERICAL verification 0.7071 = 1/rt2)
Binet's log-Gamma remainder log G(z) = (z-1/2)log z - z + (1/2)log(2 pi) + mu(z),
mu(z) = 2 int_0^oo arctan(t/z)/(e^{2 pi t}-1) dt. We PROVE |mu(z)| <= 1/(12 |z| cos(arg z))
from |arctan(t/z)| <= t/(|z| cos(arg z)) (via |z^2+s^2| >= (|z|^2+s^2) cos(arg z)) and
int_0^oo t/(e^{2 pi t}-1) dt = 1/24; on the ray arg z = pi/4 the arb computation of mu gives
ratio <= 0.7071 = 1/rt2 for x in [0.1, 1e6] (step_08/09 output). The ReE bound and the
|3/2-(u+1)log(1+3/(2u))| <= 3t/4 + (459/24) t^3/(1-3t) Taylor-remainder bound are verified
sharp (ratios 0.99..1.00, step_09 output).

## 3. [100,400] tightening (Task 2) — PROVED (arb)

29 cells, ratio <= 1.05 (`step_15.py`), each an arb enclosure from step_38 (J=24, nsub=4).
C rises smoothly from 125.69 ([100,105]) to the sup **194.3809** at [373.3456, 392.0129],
last cell [392.01,400] = 169.66. **sup C = 194.3809 < 200**, so step_34's midpoint+reach W3
port is NOT needed (conditional task, skipped as instructed). Dominant piece: W3 (~w^2 growth).
Full cell-by-cell table: stdout of step_15.py / tight_100_400.json.

## 4. Cover of [400, 10^10] (planner's sweep — used, not redone) — PROVED (arb)

43 cells, EXACTLY contiguous [400, 10^10] (verified: data[i+1][0] == data[i][1] for all i),
sup **617117.5536** at [6.6344e9, 9.9516e9]. Provenance spot-checked by re-running
`certified_y38(ys(N2), ys(N1), J=24, nsub=4)` on cells 0, 20, 41: exact match to all printed digits
(`step_16.py`). Archived copy: `sweep_cells_400_1e10.json`.

## 5. Evidence table (arb enclosures, ctx.dps=70)

| quantity | value | radius | digits |
|---|---|---|---|
| final cell C (J=192) | 49.56259362936965654550929024126328268240614812466812154794657124742081 | 2.7e-69 | 68 |
| final cell C (J=48) | 53.61568583421857309778589251677539452160470915959575395741864585441471 | 2.8e-69 | 68 |
| W1reg | 7.81474850984005936973364394174617069539632077035935395426342088504231 | 2.4e-69 | 68 |
| lamRreg | 1.390087650764052101493639038554983784331824575282770368044796567861941 | 6.5e-70 | 69 |
| sup [100,400] | 194.3808878592793 (float of arb upper) | — | — |
| sup [400,1e10] | 617117.5535603258 (float of arb upper) | — | — |
| z-factorisation ratio | 1.0000000000000000000000000000000000000000000000000000000 | 1.3e-57 | 56 |
| mu-bound worst ratio | 0.7071067811865476 (= 1/rt2 asymptotic) | — | NUMERICAL |
| g0 Binet bound ratio x=8e4 | 0.999994 | — | NUMERICAL support |

## 6. Rerun commands

  PY='<home>)/.venv/bin/python'
  cd research/08_theory/D13h_final_cell
  $PY step_38.py        # base machinery: point evals + 4 cells [100,400]         (<1 s)
  $PY step_10.py        # FINAL CELL J=48 : C = 53.6157 for all N >= 1e10         (<1 s)
  $PY step_10_J192.py   # FINAL CELL J=192: C = 49.5626 for all N >= 1e10         (~2 s)
  $PY step_15.py        # [100,400] ratio-1.05: 29 cells, sup 194.3809            (~1 s)
  $PY step_16.py        # spot-check planner sweep cells 0,20,41 (exact match)
  # planner's sweep (not redone): .../d13g_rerun/sweep_cells.json -> ./sweep_cells_400_1e10.json

Files: step_38.py (D13g base), step_08.py/step_09.py (identity checks), step_10.py + step_10_J96.py +
step_10_J192.py (final cell), step_13.py/step_14.py (W2/W3 validation vs true integrals),
step_15.py (tightening), step_16.py (sweep spot-check), step_17.py (junction + refinement),
final_cell.json, final_cell_J96.json, final_cell_J192.json, tight_100_400.json, sweep_cells_400_1e10.json.

## 7. Labels
- THEOREM (C = 617117.6, all N >= 100): **PROVED (arb)**, with the single classical citation of
  Binet's log-Gamma remainder formula (itself numerically verified to ratio 0.7071 on the ray arg z = pi/4).
- The mu-inequality |mu(z)| <= 1/(12|z|cos(arg z)), all Taylor/log remainder bounds, all monotonicity
  lemmas M1-M4: proved in the derivation; every inequality between arb enclosures machine-checked.
- NUMERICAL (support only, not used in the certificate): pointwise trend C(N): 155.5 (400) -> 83.3 (1e5)
  -> 57.5 (1e7) -> 51.2 (1e9) -> limit ~49, consistent with the final-cell 49.5626.
- Optional future tightening (NOT done, not needed): a ratio-1.05 cover of [400,10^10] (~330 cells,
  ~20 s) would lower the middle sup from 617117.6 to O(10^2-10^3); the theorem stands with C = 617117.6.
"""
open('result.md','w').write(result)
print("result.md written,", len(result), "chars")
with open('progress.log','a') as f:
    f.write("turn8: result.md written. THEOREM C=617117.6 for all N>=100 (sups 194.3809 / 617117.5536 / 49.5626).\n")
