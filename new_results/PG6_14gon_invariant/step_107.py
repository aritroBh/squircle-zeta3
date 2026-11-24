result = r"""# PG6: the regular 14-gon at the diamond end — P'''(1) resolved, and the (2,1) invariant delta(N3)

Work directory: `research/08_theory/PG6_14gon_invariant/`. All numbers below were printed by code (`step_NN.py`; `verify_pg6.py`, 22 s, **OVERALL: PASS**). Labels: PROVED / EXACT / NUMERICAL(digits) / CONJECTURE / OPEN.

---

## Task A — the P'''(1) discrepancy: RESOLVED

Two PG5 candidates disagreed by 1.05e-5: polar-perimeter gold **-5.5080732424212877** vs the exact-integrand chain 84·L·∫E₃ = **-5.5080836980534783**.

**The exact-integrand chain value is CORRECT: P'''(1) = -5.5080836980534783.** A third independent method (high-precision tanh-sinh perimeter on the *centered* chamber [-π/14, π/14], 9-point central FD + Richardson extrapolation over h = 0.02…0.0025) converges to **-5.508083698053478284**, matching the chain to ~24 digits. A roots-of-unity Fourier filter on the *full* (sliver-free) chamber agrees to ~24 digits.

**The gold value is WRONG, and the bug is located.** PG5's `step_06_gold14.py`/`step_06b_gold14.py` integrate the perimeter over the chamber `[π/2+1e-9, π/2+π/7−1e-9]`, deleting 1e-9 slivers at the two **vertex** endpoints. Near a vertex the p-ball integrand has a boundary layer: the term |cos(θ−π/2)|ᵖ with θ−π/2 ∈ [0,1e-9] multiplies each p-derivative by ≈ ln(1e-9) = −20.7, so the deleted sliver contributes an error err(p) to perim(p) with err'''(1) ≈ **+1.05e-5**. Because this error is *systematic* (the same smooth function at every p-sample), the Fourier filter at both radii (0.02 and 0.01) and the FD stencil all extract the **same** err'''(1), so they agree with each other to 7e-18 while all being wrong by 1.05e-5. Verified directly: the gold chamber (slivers deleted) reproduces −5.5080732424 to ~1e-10; the full/centered chambers give the truth −5.508083698. (NUMERICAL, 24 digits.)

**Control (12-gon, P'''(1) known = -5.2032404295176662).** The corrected pipeline (chamber *split* at the kink θ=0 — the 12-gon has a normal at π/2 so |sinθ|ᵖ vanishes at chamber center; no slivers) gives P(1)=12√(26−15√3) to 5e-70, Fourier P'''(1) to ~20 digits, FD to 2e-11. (The 14-gon has no normal at π/2, so its chamber center is smooth and needs no split.)

---

## Task B — the (2,1) invariant delta(N3): **NONZERO** (like the octagon and 12-gon; not zero like the decagon)

Ported the 12-gon engine to F = Q(2cos(2π/7)) exact arithmetic (`engine14.py`: FT field arithmetic, universal Jloc/Jtab tables, pair M1–M4/MixedW and TRIPLE recursions). Validated against quadrature to ~1e-59 on 13 integral types. The 14-gon edge data (from PG5's geom14/letters14): x = {1−α, 2−α−α², 2−α², 3−α−α²}, log(1−x) = {la, lu, −lu, −la}; Landen pairs {1↔6, 2↔5}; letters {2, α, u = α²+α−1}.

**Assembly (end-to-end validated).** E3n14 (119 term-log monomials, degrees 1–3) → 811 reduced monomials (weights {3:660, 2:134, 1:17}) → S (132 terms). **84·L·S = P'''(1) to 4.4e-50** — the whole engine+assembly is validated numerically.

**Symbol of the weight-3 non-local part.** σ(S) has 8 terms, all in (span{α,u})³ (the xᵢ, 1−xᵢ have no factor of 2). Atom symbols: Li₃ = z⊗z⊗(1−z); A(x) = Li₃(x)−Li₃(1−x)+ζ(3)+products (derived exactly, validated 1e-61); T3p = −(1/xᵢ)N with N = 2Li₃(c)−2Li₃(cX)+products (PGEN1R's repair, validated 1e-61); T3a via the S3 Goncharov recursion on numerically-verified G-words; Lam via DOD3's S_Li21 (converted −τ into PGEN1R's convention).

**Reduction (alphabet rule for the 14-gon).** σ(S) reduces to the 24 cross-ratio Li₃ symbols mod shuffle-decomposables with **EXACT residual 0** in all three F-components (rank[Li₃|Decomp] = 21; the Li₃ span mod Decomp is 2-dimensional). This is the alphabet rule confirmed at s = 5 far kink lines (14-gon). (EXACT.)

**The invariant.** δ(N₃) = Σ c_z·δ_Li₃(z), δ_Li₃(z) = (z∧(1−z))⊗z. δ **vanishes on all 22 independent Li₃-relations** (so it is well-defined / functional-equation-invariant); it is identical across 4+ particular solutions, across a minimal 2-element basis, and across both Lam-sign conventions. (EXACT, symbol level.)

\[
\boxed{\ \delta(N_3)\ =\ c_\alpha\,(\alpha\wedge u)\otimes\alpha\ +\ c_u\,(\alpha\wedge u)\otimes u\ }
\]

with (α = 2cos(2π/7), u = α²+α−1 = 2cos(π/7)):
- **c_α = (601 − 182α − 243α²)/32**  ≈ −0.118909811681569
- **c_u = (613 − 118α − 292α²)/32**  ≈ +0.369019755994578

**δ(N₃) ≠ 0** (both coefficients nonzero; wedge³ = 0 trivially). So the 14-gon's far kink lines add **genuine depth-2 (Li_{2,1}) trilogarithms** — the (2,1) obstruction is nonzero, as for the octagon and the 12-gon, in contrast to the decagon (where it vanishes). Derivation status: **EXACT / PROVED at symbol level** (the engine is exact; σ→Li₃ reduction has residual exactly 0; δ is relation-invariant; atom-reduction identities are classical/proved, numerically cross-checked to 1e-61).

**C₇ check (paper's local law).** The local ζ(3) part of S gives **84·L·z3_local = C₇ = −0.6838425376363484** (C_m = 48m·sin⁴(π/2m)cos(π/2m)(3sin²(π/2m)−1) at m = 7) to **3.6e-69**. Notably the ζ(3) coefficient in S is *entirely* local — the single/pair/triple far families contribute exactly 0 to it. (EXACT.)

---

## Task C — closed form of P'''(1): OPEN

Not attempted (δ(N₃) was the goal and is complete). The full DOD5-style closure (products + invisible constants, e.g. via the symbol-reduction c_z plus the product bookkeeping) remains open; the symbol-level reduction σ(S) = Σ c_z Li₃(z) mod products is in hand (`delta14.pkl`, `q3setup.pkl`).

---

## Evidence table (all printed by code)

| check | result | label | step/verify |
|---|---|---|---|
| centered-chamber Richardson FD P'''(1) | −5.508083698053478284 | NUMERICAL 24d | A1 / s_A2 |
| 84·L·∫E₃n (exact chain) | −5.508083698053478284 | NUMERICAL 24d | A2 / s_A4 |
| full-chamber Fourier (ρ=0.015) | agrees to ~24 digits | NUMERICAL | s_A4 |
| gold chamber (slivers) reproduces −5.5080732424 | bug located | NUMERICAL | A3 / s_A3 |
| 12-gon control P(1), Fourier P'''(1), FD | 5e-70 / ~20d / 2e-11 | NUMERICAL | s_A6 |
| engine14 vs quadrature (13 types) | ≤ 1.2e-58 | NUMERICAL | s_B5 |
| 84·L·S (assembly) vs P'''(1) | 4.4e-50 | NUMERICAL | B1 / s_B6 |
| σ(S)→24-CR Li₃ mod Decomp | residual exactly 0 (3 comps) | EXACT | B3 / s_B14 |
| δ vanishes on 22 Li₃ relations | well-defined | EXACT | B4 / s_B15,16 |
| δ(N₃) = c_α(α∧u)⊗α + c_u(α∧u)⊗u | c_α=(601−182α−243α²)/32, c_u=(613−118α−292α²)/32 ≠ 0 | EXACT | B5,B6 / s_C8 |
| δ basis/convention robustness | minimal basis = 24-CR; both Lam conventions nonzero | EXACT | s_C9,C10 |
| 84·L·ζ(3)_local = C₇ | 3.6e-69 | EXACT | B2 / s_C2 |

## Files in `research/08_theory/PG6_14gon_invariant/`
`engine14.py` (FT engine), `sym14.py` (symbol engine), `mons14.pkl` (811 monomials), `S14.pkl` (assembly), `sigma14.pkl` (weight-3 symbol), `delta14.pkl`,`delta_FINAL.pkl`, `q3setup.pkl`, `atoms_num.pkl`, `constvec14.pkl`, `geom14.pkl`/`E3n14.pkl`/`letters14.pkl` (from PG5), `Jloc.pkl`/`Jtab_t.pkl` (engine tables), **`verify_pg6.py`** (recomputes all, 22 s, OVERALL: PASS), `progress.log`, `summary.md`, `result.md`.
"""
open('result.md','w').write(result)
print("result.md written,", len(result), "chars")
open('progress.log','a').write('PG6 s_F1: result.md written; verify_pg6 PASS\n')
