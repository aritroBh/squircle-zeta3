result = r"""# E5e — non-local ζ(5) of the b = 0 hexagon arcs: PSLQ at 250–310 digits. RELATION FOUND AND CERTIFIED.

**Verdict (task branch 2): `X = P^(5)_arc/√5` lies in the ℚ-span of the 20-element weight ≤ 5 alternating-MZV
basis, with a primitive integer relation of maxcoeff = 2,574,432,000.** The relation is returned **bit-identical
by PSLQ at 250, 280, 310 (and 340) digits**, its residual at 340 digits is `5.69e-330 ≤ maxcoeff·10^(−330) =
2.57e-321` (acceptance criterion met, margin 10^8), and an independent flint LLL (δ = 0.99, 300-digit lattice)
returns **exactly ± the same vector** as its shortest row — the next-shortest lattice vector sits at the
Dirichlet scale ~5×10^14, five orders of magnitude larger. This resolves E5d's failure at 187 digits: the true
maxcoeff 2.57×10^9 was above E5d's detection ceiling ~8×10^8.

## 1. Data check (step_05, step_06)

- `../E5d_nonlocal_pslq/b0_p5_convB_350.json`: keys `e1` (r = 2), `e2` (r = 1/2), 337 decimal digits each.
- **e1, e2 agree through all 336 first decimals** (common prefix 341 of 342 chars; last stored digit differs:
  `…5285` vs `…5289`); `e1 − e2 = 4.0e-337` absolute = `4.69e-340` relative — consistent with E5d's convergence
  log (edge 1: 9.49e-339; edge 2: 6.43e-339). ⇒ ~336 reliable digits; prefix matches E5d's 187-digit value.
- `X = P^(5)_arc/√5 = −381.140137469830729527176573357265848015393511891093057029…`
- Basis recomputed at 360 dps from `../E5c_nonlocal_z5/mzv.py` (20 elements; stored in `consts360.json`):
  `1, log2, π², log2², ζ3, π²log2, log2³, π⁴, ζ3log2, π²log2², log2⁴, Li4(1/2), ζ5, ζ2ζ3, π⁴log2, ζ3log2²,
  π²log2³, log2⁵, Li5(1/2), Li4(1/2)log2`.

## 2. PSLQ runs (step_08): [X, basis] ∈ ℝ^21, tol = 10^(−(d−25)), maxcoeff ≤ 1e14

| precision | tol | result | residual at working precision | maxcoeff |
|---|---|---|---|---|
| 250 digits | 1e-225 | rel R (below) | −8.7873e-244 | 2,574,432,000 |
| 280 digits | 1e-255 | **same R** | −1.6175e-272 | 2,574,432,000 |
| 310 digits | 1e-285 | **same R** | +2.5154e-302 | 2,574,432,000 |
| 340 digits (step_10) | 1e-315 | **same R** | +5.6841e-330 | 2,574,432,000 |

`R = [7680000, 0, 1175040000, 49152000, 2574432000, −501811200, −111744000, 927640800, 6229440, 397036800,
45175200, 287114250, 381542400, 77518080, 52715520, −1240896, −1751040, −10005440, −2977497, −24330240,
−56401920]` (coefficients of `[X, 1, log2, …, Li4(1/2)log2]`); gcd = 1 (primitive).

Residuals scale like `maxcoeff·10^(−dps)` (≈2.6e-241 / 2.6e-271 / 2.6e-301 predicted) — the true-relation
signature, and ~10^19 below tolerance (contrast E5d's spurious tolerance-level relations).

## 3. Certification at 340 digits (step_09, step_10)

| check | value | bound | pass |
|---|---|---|---|
| residual, e1 data, mpmath 340 dps | 5.6841e-330 | maxcoeff·1e-330 = 2.574432e-321 | ✓ (margin 4.5e8) |
| residual, **e2 leg** (independent evaluation) | 4.2853e-330 | same | ✓ |
| residual, **exact rational** arithmetic on the 350-digit truncated decimals | 5.6886e-330 | same | ✓ |
| PSLQ relation identical at 250/280/310/340 digits | True | required | ✓ |

## 4. LLL cross-check (step_11, flint `fmpz_mat.lll`, δ = 0.99)

Lattice rows `(e_i | round(v_i·10^300))` ∈ ℤ^22, v from the 350-digit data. `‖R‖ = 3.088e9` vs Dirichlet scale
10^(300/21) ≈ 1.93e14. LLL row 0 = **−R exactly** (norm 3.136e9, last coordinate 5.46e8 ≈ rounding noise);
rows 1–5 have norms 5.29e14 … 6.62e14. The relation is the unique short vector, separated by 5 orders of
magnitude.

## 5. Exact closed form (task item 2)

`P^(5)_arc = √5·X` with

```
X = −153·log2 − (32/5)·π² − (26817/80)·log2² + (3267/50)·ζ(3) + (291/20)·π²log2 − (386517/3200)·log2³
    − (6489/8000)·π⁴ − (20679/400)·ζ(3)log2 − (18823/3200)·π²log2² − (382819/10240)·log2⁴ − (1242/25)·Li4(1/2)
    − (20187/2000)·ζ(5) − (858/125)·ζ(2)ζ(3) + (6463/40000)·π⁴log2 + (57/250)·ζ(3)log2²
    + (31267/24000)·π²log2³ + (992499/2560000)·log2⁵ + (396/125)·Li5(1/2) + (918/125)·Li4(1/2)log2
```

(equivalently `−7680000·X = Σ R[i]·b_i`; denominator 7680000 = 2^12·3·5^4.)

**340-digit re-check (step_12):** closed form vs data:
`X − closed = 7.37e-337`; `P5 − √5·closed = 1.66e-336` — agreement to the data's own ~336-digit accuracy. ✓

## 6. ζ(5) and ζ(2)ζ(3) coordinates vs the local law (step_12)

In THIS basis (coords of `X`; coords of `P^(5)_arc` are √5 × these):

| constant | coord in X | local-law coord | non-local part |
|---|---|---|---|
| ζ(5) | −20187/2000 = −10.0935 | −2079/250 = −8.316 | N^(5) = −20187/2000 + 2079/250 = **−711/400**, i.e. `N^(5) = −711√5/400 ≈ −3.9746108300058761854` for P |
| ζ(2)ζ(3) | −858/125 = −6.864 | −477/250 = −3.816 | N^(2,3) = −858/125 + 477/250 = **−1239/250**, i.e. `N^(2,3) = −1239√5/250 ≈ −11.081952896488957735` for P |

Both non-local parts are **nonzero and large**: the arc value is far from the local law, with rational (in
X-space) defect −711/400 on ζ(5) and −1239/250 on ζ(2)ζ(3).

**Basis-dependence:** the coordinates of `Li5(1/2)` (+396/125), `Li4(1/2)` (−1242/25) and `Li4(1/2)log2`
(+918/125) are nonzero. Changing the weight-5 basis (e.g. eliminating Li5(1/2) via the standard
Li5(1/2)-identity, or swapping to ζ(4,1)-type elements) reshuffles the ζ(5) and ζ(2)ζ(3) coordinates, so
N^(5), N^(2,3) are basis-dependent numbers, not invariants — as the task states.

## 7. Consistency with E5d

E5d (187 digits) proved any relation has maxcoeff ≳ 2×10^8 and could detect nothing above ~8×10^8; the true
maxcoeff is 2.57×10^9, exactly in the predicted gap, and the E5d estimate "Detecting maxcoeff ~1e12 needs
~250–270 digits" is borne out (relation already locked at 250 digits; needed only 2.6e9).

## 8. Files (work dir `…/08_theory/E5e_pslq_340/`)

- `step_05.py` — json data inspection (e1/e2 agreement, 336 digits)
- `step_07.py` — 360-dps constants & basis → `consts360.json`
- `step_08.py` — PSLQ at 250/280/310 digits → `pslq_runs.json` (identical R)
- `step_09.py` — 340-dps residual (e1+e2), gcd, exact rational coordinates → `relation.json`
- `step_10.py` — exact-rational residual; PSLQ at 340 digits (identical R)
- `step_11.py` — flint LLL cross-check (row 0 = −R exactly)
- `step_12.py` — closed-form re-check at 340 digits; N^(5), N^(2,3)
"""
open("result.md","w").write(result)
print("result.md written,", len(result), "chars")
print(open("result.md").read()[:600])
