DOD5 complete — termwise closure of all 36 depth‑2 constants, the corrected assembly, and the closed form of the 12‑gon's P‴(1), verified to 611 digits.

# DOD5: the regular 12‑gon's P‴(1) in closed form

## 1. What was wrong in DOD4, and the two bugs found here

**(a) The one‑shot Ξ36 PSLQ failure is a coefficient‑size explosion, not a convention or pairing error.** The atom/coefficient pairing is correct: Σc^Lam·Lam + Σc^T3a(Swγ−Swδ)/xᵢ with DOD2's exact ℚ(√3) coefficients equals Ξ − Σc·R to **1.37e‑614** (printed). With the corrected prefix symbol, each of the 36 constants closes separately (below). Summing the 36 individual remainders with the 1‑part coefficients reproduces DOD4's rem0 = −579.19… (I get −579.1928… vs DOD4's −579.1914…, the small gap being freedom in the product decomposition) — but the combined relation's rational coefficients involve **94 distinct denominators** from the 36 relations, hopelessly beyond PSLQ at maxcoeff 10¹³. Termwise closure is the fix.

**(b) Correction to DOD2/DOD3: the "~1e‑4 branch terms" in 3–4 R pieces were quadrature error, not branch terms.** DOD2's `gdswaps620.pkl` computed the 8 swaps with a γ or δ argument > 1 by plain `mp.quad` over an interior log singularity; those 8 values are wrong at up to **9.02e‑5**. Recomputed by tanh‑sinh quadrature with explicit splitting (stable: |620−680 dps| = 1.2e‑618, |680−740| = 4e‑678). Consequence: **all 12 symbolic R pieces (DOD3's `build_R`) match T3a − (Swγ−Swδ)/xᵢ to ≤ 1.2e‑619** — no branch‑term products are needed anywhere. Also, where the prefix series converges with y > 1, it agrees with the real‑analytic value to ~7e‑622 (no branch ambiguity).

## 2. Task 1 — the 36 constants (all printed with exact arguments)

Convention: prefix f(x,y) = Σ_{0<k<m} x^m y^k/(m²k) = Goncharov's Li_{1,2}(y,x). 12 Lam(i,j) = f(xᵢ, xⱼ/xᵢ); 12 SwG(i;j,k) = f(γᵢₖ/γᵢⱼ, γᵢⱼ); 12 SwD(i;j,k) = f(δᵢₖ/δᵢⱼ, δᵢⱼ), γᵢⱼ = xⱼ/(xⱼ−xᵢ), δᵢⱼ = γᵢⱼ(1−xᵢ), i,j,k ∈ {1,2,5,6}, j<k. Full exact table printed (e.g. Lam(1,2): x=(3−2√3)/3, y=3+√3, value 0.0231418623633688955…; SwG(6;1,5): x=1+√3, y=4−2√3, value −1.16512707002030991…). Values: series where |x|<1∧|xy|<1∧|y|≤1 (14), else real‑analytic C_real with tanh‑sinh (22). All 36 verified against series/quadrature at well over 100 digits.

## 3. Task 2 — all 36 constants close separately (PROVED symbol + NUMERICAL remainder)

Symbol S_prefix(x,y) = S(Li₂Li₁) − S_Li21(x,y) − S_Li3(xy) (DOD4's, from the stuffle). Exact decomposition over the 8‑arg Li₃ basis + 19 products (tensor residual exactly 0, asserted in code — **PROVED**); remainder closed in invisible10 = {X, √3X : X ∈ ζ(3), L(3,χ₁₂), ζ(2)log2, ζ(2)log3, ζ(2)log(2+√3)} by PSLQ at 620 **and** 450 dps, accepting only identical relations. **All 36 pass**: max coefficient 3024, residuals ≤ 7.7e‑618 at the 620‑digit working precision (**NUMERICAL**, two precisions). Saved in `dod5_closures36.pkl`.

## 4. Tasks 3–4 — assembly (all printed)

- **Ξ36** = Σc·(36 closures) matches the 36‑term value 0.01383198578110140502410683742522418757579… to **9.17e‑619** (8 Li₃ + 19 products + invisible, exact ℚ(√3) coefficients; Li₃ coefficients equal DOD4's split values, confirming the symbol side was right).
- **R pieces**: Σc^T3a·R_sym rebuilt exactly, identical to DOD3's `R_li3_content.pkl` (12 Li₃ points); all 12 verified numerically (≤1.2e‑619).
- **Ξ = Ξ36 + ΣcR** = −0.01440930466893363902462… matches DOD2's `Xi_from_S` to **1.37e‑614**.
- **Depth‑1 part**: exact ℚ(√3)(+log) coefficients of all 60 atoms extracted (`dod5_atom_coefs_all.pkl`; Lam/T3a match DOD2). S_d1 = Σc·atom + singles matches S_nod2 to **5.4e‑614**. A‑atom Nielsen formula verified 1.5e‑620; my simpler T3p closed form (1/p)[A²log u − 2A·Li₂(γu) + F(γu)]_{1−p}¹ verified vs all 12 atoms to 3.5e‑620.
- **P‴(1) = 72·L·(S_d1 + Ξ)** = −5.20324042951766621029582970427243461666280533… vs `P3_630.pkl`: **|diff| = 2.59e‑612 (611.6 digits)** — well beyond the ≥500‑digit requirement.

## 5. Task 5 — merged form, L(3,χ₁₂), ζ(3)

23 raw Li₃ points in S; all reduce over the 8‑arg basis (symbol‑exact + two‑precision remainders, all pass, ≤1.2e‑619). Full bookkeeping closes: S_calc = S_full to **1.37e‑614**.

**Final merged form:** S = Σ₈ Tot3ⱼ·Li₃(zⱼ) + Σ₁₉ products + Σ₅ invisible + (named depth‑1 product atoms), with **all 8 Li₃ arguments surviving** (minimum; the rank‑8 quotient is DOD3's PROVED span theorem): Li₃ at {2−√3, (√3−1)/2, √3/3, 7−4√3, (2−√3)/2, (2√3−3)/3, (3−√3)/6, 1/3} with S‑coefficients (−73/24+37√3/8), (11/2−39√3/4), (97/24−63√3/8), (19/192−13√3/64), (−1/8+9√3/8), (2/3−3√3/4), (7/4−13√3/4), (−17/96+13√3/16).

**L(3,χ₁₂) occurs**: total coefficient in S = **459/16 − 93√3/16 ≠ 0**, i.e. in P‴(1): 72L(459/16 − 93√3/16)·L(3,χ₁₂) ≈ 185.9472649216808190731434·L(3,χ₁₂) — the first Dirichlet L‑value at weight 3 in a polygon's P‴(1) (well‑defined modulo the standard independence of the 10 invisible constants; the Li₃ coefficients are symbol‑proved). Also L(2,χ₁₂) = π²√3/18 confirmed (0.0 at 120 dps), so √3·L(2,χ₁₂)=ζ(2) collapses weight‑2.

**ζ(3) vs the local law**: total ζ(3) coefficient in S = 2131/2304 − 481√3/768 (in P‴(1): −1.596572216112899742201299). The singles part carries exactly (2−3√3)/32 = 1/16 − 3√3/32, and 72L·(2−3√3)/32 = (9/4)(2−3√3)L = **C₆ = −0.9974443605274988364594569 exactly**, matching the paper's local law C_m = 48m·sin⁴(π/2m)cos(π/2m)(3sin²(π/2m)−1) at m = 6. The remainder (1987/2304 − 409√3/768) is the non‑local (basis‑dependent) part from A‑atom Nielsen terms and Li₃ remainders.

## 6. Evidence table (all printed by code)

| check | residual | digits |
|---|---|---|
| 36 constants fresh vs stored atoms/gdswaps | ≤9e‑5 → **8 gdswaps values wrong** (corrected) | 620 |
| 36 closures: symbol decomposition | tensor residual exact 0 (PROVED) | — |
| 36 closures: remainder PSLQ identical @620/@450 | ≤ 7.7e‑618, maxcoef 3024 | 620 |
| series vs real‑analytic (y>1 convergent) | 7.1e‑622 | 620 |
| 12 R sym vs accurate values | ≤ 1.2e‑619 (no branch terms) | 620 |
| Ξ36 closed vs value | 9.17e‑619 | 620 |
| Ξ closed vs Xi_from_S | 1.37e‑614 | 620 |
| S_d1 assembly vs S_nod2 | 5.41e‑614 | 620 |
| T3p formula (12), A‑Nielsen (4) vs atoms | 3.5e‑620 / 1.5e‑620 | 620 |
| 23 Li₃ reductions (two precisions) | ≤ 1.2e‑619 | 620 |
| S_calc (fully merged) vs S_full | 1.37e‑614 | 620 |
| **72·L·S_calc vs P3_630** | **2.59e‑612 (611.6 digits)** | **620** |

## 7. Files (in `DOD5_dodecagon_termwise/`)

`dod5_task1.pkl` (36 constants+values), `dod5_closures36.pkl`, `dod5_xi36_closed.pkl`, `dod5_R_content.pkl`, `dod5_all_li3_reductions.pkl` (23 pts), `dod5_tallies.pkl` (merged Tot3/TotP/TotI), `dod5_rawli3_coeffs.pkl`, `dod5_atom_coefs_all.pkl`, `dod5_xi_p3_closed.pkl`; libraries `dod3_symbol.py`, `dod3_R.py`, `dod4_lib.py` (copied); computations in auto‑saved `step_04.py`–`step_64.py` (key ones: step_28/30 (values+quadrature fix), step_32 (R‑piece resolution), step_37 (36 closures), step_41 (Ξ36), step_46 (Ξ, P‴(1)), step_50–54 (depth‑1, T3p, reductions), step_57–59 (merged tallies + global check), step_61 (Task‑1 table + L(2,χ₁₂))). **`verify_dod5.py`** (9.3 KB) is written to disk: it recomputes everything from scratch at 330 digits and re‑checks the 36 identities, 23 reductions, 12 R pieces, Ξ36, and 72L·S_calc vs P3_630 (target < 1e‑290); every constituent check in it was printed above at 620 digits, but its final end‑to‑end run was cut by the turn limit (a missing `import os` in my launcher, not in the script) — run `python verify_dod5.py` to reproduce; honestly labelled **NUMERICAL** like all remainder closures.

**Labels:** PROVED — S‑unit letters, all symbol decompositions (tensor residual 0), R symbolic pieces, exact coefficients everywhere, rank‑8 minimality. NUMERICAL — every remainder/relation closure (identical at two precisions, residual at working precision) and all value checks.