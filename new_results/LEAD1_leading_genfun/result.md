# LEAD1 — the generating function of the leading coefficients c(w) of the weight parts of M_k

## VERDICT

The leading coefficients `c(w) = [k^{w-2}] W_w(k)` for **w = 2..10** were computed from the exact top kernel coefficients `D_{w-1}` by quadrature (405 digits for w ≤ 9, 250 for w = 10). All four controls pass to ~410 digits. All of c(5)..c(9) are identified in the Gamma-world monomial basis with integer relations **identical at two precisions**. The generating function is

**c(w) = eps_{w-2} · A_{w-1} / (2 (w-2)!) · [y^w] E(y)**,

```
E(y) = Gamma(1+y/2) Gamma(1/2-y/2) / sqrt(pi)
     = 2^y Gamma(1-y) Gamma(1+y/2) / Gamma(1-y/2)              [duplication formula]
     = exp( y log 2 + sum_{n>=2} s_n y^n/n ),   s_n = zeta(n) (n even),  s_n = (1-2^{1-n}) zeta(n) (n odd),
```

with `eps_m = +1,-1,-1,+1` (period 4, `= (-1)^floor((m+1)/2)`) and `A_n` the Euler zigzag numbers (`sec x + tan x = sum A_n x^n/n!`: 1,1,1,2,5,16,61,272,1385,7936,...). Verified against **every** computed c(w): w = 2..9 to 406 digits, w = 10 to 237 digits.

Equivalent forms (since `eps_m A_{m+1}/m! = [z^m](1-sinh z)/cosh^2 z`, because `(1-sinh z)/cosh^2 z = d/dz(sech z + tanh z)`; symbolic, step_54):

- `F(u) = sum_{w>=2} c(w) u^w = (1/2) [ u^2 (1-sinh u)/cosh^2 u  ⊙  E(u) ]` (Hadamard = coefficientwise product);
- reduced GF is a pure Gamma product: `sum_{w>=2} [2(w-2)!/(eps_{w-2} A_{w-1})] c(w) u^w = E(u) - 1 - u log 2`.

## 1. Values (60 digits shown; 405/250 digits in `leading_coeffs.json`)

| w | c(w) | sign |
|---|---|---|
| 2 | 0.531346770191606965284879423243172540237375713200336006150606 | + |
| 3 | -0.463054519888073265020183366821329463829457953996109311061807 | - |
| 4 | -0.512151959461527680222200558606808431407661025164465171146113 | - |
| 5 | 0.40948119421256435691521921506651565742131739799942432899534 | + |
| 6 | 0.33568755822402035700767227984351182721153530440389179214486 | + |
| 7 | -0.253122416960440487223947539490993007431297696004986069718015 | - |
| 8 | -0.189242541413378208209141417583276128090266362460591787429044 | - |
| 9 | 0.137263340792142938015670190422387507617099134460586704785821 | + |
| 10 | 0.0984598871418917972706502940864227127442040035844511774009152 | + |

## 2. Controls (dps 410, step_14; two independent quadrature runs, steps 13/26/27)

| control | residual |
|---|---|
| c(2) − (π² + 6 λ²)/24 | 0.0 (exact to 410 digits) |
| c(3) + (ζ(3)/8 + π²λ/24 + λ³/12) | 6.2e-412 |
| c(4) + (π⁴/320 + π²λ²/48 + λ⁴/48 + ζ(3)λ/8) | 1.2e-411 |
| c(6) − (identified form /9!) | −2.5e-411 |
| run(dps 250) − run(dps 420), w=2..9 | ≤ 6.5e-238 |

## 3. Identifications (λ = log 2, z_k = ζ(k))

LLL (python-flint, `gram='exact'`) relations — **identical at dps 420 (D=390) and dps 250 (D=225)**, and equal to the law's closed-form relations (step_25):

- `3456 c(5) = 12λ⁵ + 20π²λ³ + 180λ²z3 + 9π⁴λ + 30π²z3 + 270z5` (matches w6's c₃⁽⁵⁾)
- `362880 c(6) = 168λ⁶ + 420π²λ⁴ + 5040λ³z3 + 378π⁴λ² + 2520π²λz3 + 22680λz5 + 61π⁶ + 3780z3²` (matches identify_w6.md)
- c(7): 12 monomials, max coeff 2.9e7, relation in `lll_relations_dps400.json`
- c(8): 14 monomials (incl. `21437π⁸`, 21437 = 17·1261), max coeff 2.6e9
- c(9): 19 monomials (incl. `23545 z9`, 23545 = 5·17·277), max coeff 5.9e10
- c(10) = `(31/315) e_10`, closed form in `c10_closed_form.json` (24 monomials; every coeff divisible by 31 = A₉/256 except 1891=31·61, 39091=31·1261, 527=17·31, 2449=31·79)

LLL residuals (dps 400): −1.5e-402, −6.4e-401, −6.8e-400, −5.2e-397, −2.7e-395; (dps 250): −1.5e-238, 0.0, 5.0e-234, 3.2e-230, −3.8e-227.

## 4. The exact top-coefficient table D_d (the combinatorial heart)

Extracted from `kern_k.pkl` (WTS/theory) by Newton differences over k = d+1..18; with `δ_d = [x^d](sech x + tanh x) (= eps_{d-1}A_d/d!)`:

- pole −a (all keys): `D_d = (-1)^{j+1} binom(d,j) δ_d`
- pole a (odd d only): `D_d = (-1)^j binom(d,j) δ_d` (×1/2 on the last key)
- pole −b (even d, j < d/2): `D_d = (-1)^j binom(d,j) δ_d`
- pole 0 (diagonals only): `D_d = (-1)^j binom(d,j) δ_d /2` (odd-d last key also carries −b)
- everything else: 0 (all sqrt2 parts vanish)

**EXACT** for d = 1..8 (all keys, all poles, step_36) and d = 9 (step_57; determined fit on k=10..18, zero spare points, still matches everywhere). For d ≥ 10 this is the unproved tower law T of WTS/theory.

## 5. Derivation / proof status (full account in `genfun.md`)

With `u = (v−1/v)/2` one has exactly `t = (1+u)/2`, and the pole-combo measures become (exact quotient-ring algebra, step_38): odd `[1+√2 u/S]/(1−u²)`, even `[1−√2/S]/(1−u²)`, diag `(1−u)/2 × even`, oddlast `odd − 1/(2(1−u))`, `S=√(1+u²)`.

- **Even d (odd weights): PROVED for all d given the D-table.** Parity collapse to a single convergent integral, then `∫_0^1 (−log w)^d [(1+w)/√(1+w²) −1] dw/w = d! e_{d+1}` proved for all even d via the generating function `F(y) = Σ (−1)^n binom(2n,n)/4ⁿ [1/(2n−y)+1/(2n+1−y)]` and `(F(y)+F(−y))/2 = (E(y)−E(−y))/(2y)` (Beta substitution + IBP; |y|<1). ⇒ c(3), c(5), c(7), c(9) proved given the table (c(2), c(4) were already proved independently).
- **Odd d (even weights):** the single integral `(1/(2(d+1))) ∫_0^1 (log w)^{d+1} (1+w)/(1+w²)^{3/2} dw = (d!/2) e_{d+1}` is **proved** for all odd d (same machinery). The key-sum equals this integral **checked** to 37–80 digits for d = 1..9 with real kernel data; the all-d algebraic collapse is **not proved** (finite-part obstruction in `G_odd = B(x,y)/2 + √2 2^{−(x+y)}J(x,y)`; the decomposition `(piece1 − boundary) − piece2 = (d!/2)e_{d+1}` with Bose-Einstein-log integrals `A_{p,q}` is verified to ~55 digits for d ≤ 9).

## 6. Evidence table (all numbers printed by code)

| check | step | residual / agreement |
|---|---|---|
| law vs c(w), w=2..9 (dps 410) | step_21 | ≤ 4.8e-406 |
| law rerun vs c(w), w=2..9 (dps 240) | verify_lead1.py D | ≤ 9e-228 |
| law vs c(10) from exact D₉ (dps 250) | step_60 | −4.6e-237 |
| E(y) series == Γ product (y=0.1, 13 terms) | verify D0 | −1.1e-31 (truncation) |
| even-d single integral == (d!/2)e_{d+1}, d=0..16 | step_34, verify E | ≤ 1e-61 |
| odd-d single integral == (d!/2)e_{d+1}, d=1..11 | step_41 | ≤ 2.5e-59 |
| key-sum == single integral, d=1..8 (real data) | step_50 | ≤ 5.3e-43 |
| odd-d collapse d=9 (real D₉ data) | step_61 | 1.6e-37 (dps 80) |
| (F(y)+F(−y))/2 == (E(y)−E(−y))/(2y), y=0.3 | verify E | −4.3e-46 |
| zigzag identities + duplication formula (symbolic) | step_54 | True / 0 |

## 7. Deliverables (all in `new_results/LEAD1_leading_genfun/`)

- `result.md` — this verdict (self-written)
- `leading_coeffs.json` — values (405/250 digits), closed forms, LLL relations, law
- `genfun.md` — the closed form and its derivation, proved-vs-checked per step
- `verify_lead1.py` — reruns A (exact D-table), B (quadrature + controls), C (LLL identifications vs stored relations), D (law), E (proof-side integrals); **ALL PASS in ~20 s** (24/24 checks, residuals printed)
- Data: `Dtop.json`, `Dtop_d9.json`, `Nints_dps420.json`, `Nints_dps250.json`, `cw_dps420.json`, `cw_dps250.json`, `cw_closed_forms.json`, `c10_closed_form.json`, `lll_relations_dps400.json`
- Producing steps: step_08 (D-extract), step_13/26 (quadrature 420/250), step_14 (controls), step_21 (law 406-digit), step_22 (closed forms), step_24/27 (LLL two precisions), step_33/34/41 (integration identities), step_36 (exact D-table), step_38 (u-measures), step_47/50/56 (odd-d decomposition), step_54 (symbolic identities), step_57/60/61 (D₉ and c(10)), step_67 (final deliverables)

**Open:** the D-table for d ≥ 10 (tower law T); the all-d algebraic proof of the odd-d (even-weight) key-sum collapse (verified d ≤ 9). Everything else stated as "verified/identified/exact" above was printed by the scripts listed.