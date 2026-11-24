# DOD2 — finishing the regular 12-gon: does depth 2 survive in P'''(1)?

**Work directory:** `new_results/DOD2_dodecagon_finish/` — `result.md` (this report), `step_01.py`–`step_17.py`, pickles listed below. All numbers below were printed by code.

## Bottom line (honest verdict, Task 4)

**Depth 2 survives at m = 6 (numerical, coefficient bound ~6·10⁶ at 620 digits).** The exact depth‑2 remainder Ξ of the 12‑gon's S (validated to 1.4e‑614 against the full symbolic assembly) does **not** reduce to the full 102‑element depth‑1 basis — its 1‑part Ξ_p and √3‑part Ξ_q each show **no genuine integer relation with maxcoeff ≲ 6·10⁶ at 620 digits** (the effective limit for a 103‑dimensional lattice at this precision, ~10^(620/103)). Individual depth‑2 constants were already excluded to 10¹³ by DOD1. The machinery was validated by two positive controls recovered exactly (CTL relation, maxc 2; R_p relation, maxc 6). Caveat, stated plainly: a reduction with coefficients **larger** than ~6·10⁶ cannot be excluded at 620 digits; the planned decisive 1350‑digit run (bound ~10¹³) and PRV1 derivative‑method attempt were not completed within the turn budget. Depth ≥ 3 remains excluded structurally (DOD1).

## Task 1 — the exact depth‑2 remainder Ξ (COMPLETE)

From the validated assembly (S is linear in all named atoms; 72·L·S = gold P3_630 to 8.95e‑611), the depth‑2 remainder is

**Ξ = Σ_{12 pairs} c^Lam_{ij}·Lam(i,j) + Σ_{12 triples} c^T3a_{i;jk}·T3a(i;j,k)**, Lam(i,j) = Li₂,₁(xᵢ, xⱼ/xᵢ),

with exact ℚ(√3) coefficients (extracted from S_atom_decomp, `radsimp`, verified 0–8e‑120 at 120 dps):
- Lam: (1,2),(1,5),(5,2) = (9−3√3)/32; (1,6),(5,1),(5,6) = 3√3/32; (2,1),(2,6) = 3/32; (2,5) = (−3+3√3)/32; (6,1) = 3/16; (6,2),(6,5) = (−3+3√3)/16.
- T3a (i;jk): (1;25) = −1047/16 + 1209√3/32; (1;26) = (4827−2787√3)/64; (1;56) = (−8349+4821√3)/64; (2;15) = −2649/32 + 765√3/16; (2;16) = (6117−3531√3)/64; (2;56) = 3057/32 − 441√3/8; (5;12) = −1527/32 + 441√3/16; (5;16) = (−6117+3531√3)/64; (5;26) = 1767/32 − 255√3/8; (6;12) = (4173−2409√3)/64; (6;15) = (−7239+4179√3)/64; (6;25) = 1047/16 − 1209√3/32.

**Numerically Ξ = −0.014409304668933639025…, Ξ_p = −1.7708150073548233046, Ξ_q = 1.0140613052518921951** (Ξ = Ξ_p + √3·Ξ_q); |Ξ(S−S_d1) − Ξ(extracted)| = **1.365e‑614** at 620 digits (step_11).

**Equivalent 36‑constant form.** With the corrected T3a reduction (below), each T3a(i;j,k) = (1/xᵢ)[Swpγ − Swpδ] + (depth‑1), Swpγ = Li₂,₁(γᵢₖ/γᵢⱼ, γᵢⱼ), Swpδ = Li₂,₁(δᵢₖ/δᵢⱼ, δᵢⱼ), so Ξ ≡ 12 Lams + 24 γ/δ‑cross swaps (mod depth‑1) — values in `gdswaps620.pkl`; the swap formula SwpR(a,b) = C_real(a,b) + ln|1−a|·(Li₂R(b) − Li₂R(b/a)) − Li₃R(b) was verified against DOD1's 12 pair‑swaps (worst **0.0**) and against the Li₂,₁ series in 8 convergent γ/δ cases (worst **2.361e‑621**).

**Correction to DOD1 (found here).** DOD1's step_71 Q(a,b) antiderivative is missing a term: the primitive of the A′B′ constant piece is **A·B·ln u**, not A·ln u (error ~0.2845 on (2,5,6) at 150 dps). After the fix, and after replacing the complex Cdef by the real C_real with singularity splitting, the closed form T3a = [−ln(1−xᵢ)³/3 + M(a) + M(b) + Q(a,b)]/xᵢ verifies against quadrature on **all 12 triples, worst 5.232e‑150** (150 dps), and exactly (0.0) at 620 digits on (2,5,6). Consequently each R(i;j,k) = T3a − (Swpγ−Swpδ)/xᵢ is **provably depth‑1**: the fully tracked symbolic expansion (all logs factored exactly into {l2,l3,lu} over ℚ(√3); Li₂/Li₃ at γ/δ, cross‑ratio and 1/γ points, reducible by Landen/reflection/inversion) matches R_quad to **1.926e‑621**, and R(2,5,6)/(1+√3) reduces with maxc 6, residual **9.59e‑622**.

## Task 2 — Ξ does not reduce to depth 1 (at bound ~6·10⁶)

Basis (102 elements, as briefed): ζ(3); ζ(2)·{l2,l3,lu}; 10 log cubes; Li₃(x₁,x₂,√3/2,√3/3); Li₂(x₁),Li₂(x₂)×3 logs; Li₃ at 12 γ/δ points; Li₂ at 6 γ/δ reps ×3 logs; Li₃ at 12 (1−1/γ) points; Li₂ at 12 (1−1/γ) points ×3 logs. Internal relations quotiented (22 found; elim reconstruction 2.85e‑618). LLL (sympy exact `Matrix.lll`, C = 10⁵⁴⁰, 620 dps) on [B₁₀₂ | Ξ_p] and [B₁₀₂ | Ξ_q]: 81 target‑involving vectors each, **all Diophantine noise** (maxc ~5–9·10⁶, residuals ~1e‑535 vs true‑relation floor maxc·1e‑615) → **no genuine relation with maxcoeff ≲ 6·10⁶**. Positive controls in identical lattices: CTL relation found (maxc 2, resid 1.18e‑621) and R_p found (maxc 6) — so the test finds relations when they exist.

## Task 3 — status

Ξ does not reduce at the achieved bound, so no depth‑1 closed form for P'''(1) exists at this precision level; the surviving directions are exactly Ξ_p, Ξ_q (equivalently the 36‑constant swap/Lam combination). What remains for full settlement: (a) recompute atoms/basis at ~1350 digits (timed: 23 s/T3a‑quad, ~40 s/C_real‑quad with splitting — feasible) and rerun the LLL to push the bound to ~10¹³; (b) PRV1's derivative method on the surviving combination. Not done (turn limit).

## Evidence table (all printed by code)

| check | residual | digits |
|---|---|---|
| 24 ℚ(√3) atom coefficients vs original | ≤ 7.8e‑120 | 120 |
| Ξ(extracted) vs Ξ(S−S_d1) | 1.365e‑614 | 620 |
| Ξ − (Ξ_p + √3 Ξ_q) | 9.207e‑621 | 620 |
| swap formula vs DOD1 pair‑swaps (12) / vs series (8) | 0.0 / 2.361e‑621 | 620 |
| DOD1 step_71 Qab bug (unfixed (2,5,6)) | 0.2845 | 150 |
| corrected T3a_cf vs quad, all 12 triples | 5.232e‑150 | 150 |
| corrected T3a_cf vs quad, (2,5,6) | 0.0 | 620 |
| R(2,5,6) tracked explicit vs quad | 1.926e‑621 | 620 |
| R_p = R(2,5,6)/(1+√3) reduces | maxc 6, resid 9.59e‑622 | 620 |
| CTL control relation | maxc 2, resid 1.18e‑621 | 620 |
| basis quotient (22 rels) reconstruction | 2.85e‑618 | 620 |
| S_full vs stored S_pq (= gold) | 8.953e‑611 | 620 |
| **Ξ_p, Ξ_q vs 102‑basis: NO genuine relation, maxcoeff bound** | **~6·10⁶** | **620** |

## Files

`step_01` (explore), `step_02`/`step_08` (coefficient extraction, `atom_coefs_pq.pkl`), `step_03` (24 γ/δ swaps, checks, `gdswaps620.pkl`), `step_05`/`step_06` (102‑basis + quotient, `basis102_620.pkl`, `basis_quotient.pkl`), `step_07` (first Ξ test, `Xi_620.pkl`), `step_09`/`step_12`/`step_13` (diagnostics), `step_10` (CTL control), `step_11` (Ξ validation, `Xi_verified_620.pkl`), `step_14` (T3a_cf fix at 620), `step_15` (exact log factorization), `step_16` (tracked R expansion, `R256_sym.pkl`), `step_17` (Ξ verdict LLL, `xi_test620.pkl`). Reused from DOD1: `dod_engine.py`, `P3_630.pkl`, `atoms620.pkl`, `S_assembly.pkl`, `S_atom_decomp.pkl`, `S_pq_620.pkl`, `gamdel620.pkl`.

**Verdict for the paper:** *depth 2 survives at m = 6* — numerical, excluded coefficient bound ~6·10⁶ at 620 digits for the combination (10¹³ for the individual constants); not a proof of irreducibility, and a higher‑precision (≥1350‑digit) run is the identified next step to make the bound decisive.