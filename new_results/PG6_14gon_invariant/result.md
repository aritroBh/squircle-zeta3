# PG6: the regular 14-gon at the diamond end — P‴(1) resolved, and the (2,1) invariant δ(N₃)

All numbers below were printed by my own code. `verify_pg6.py` recomputes everything in **22 s → OVERALL: PASS**. Labels: PROVED / EXACT / NUMERICAL(digits) / OPEN.

---

## Task A — the P‴(1) discrepancy: RESOLVED

**The exact-integrand chain value is CORRECT:  P‴(1) = −5.5080836980534783.** The polar-perimeter gold (−5.5080732424212877) is WRONG, and the bug is located.

- A **third independent method** — high-precision tanh-sinh perimeter on the centered chamber [−π/14, π/14], 9-point central FD + Richardson extrapolation (h = 0.02…0.0025) — converges to **−5.508083698053478284**, matching the chain to ~24 digits. A sliver-free roots-of-unity Fourier filter (ρ=0.015) agrees to ~24 digits.
- **The bug:** PG5's `step_06_gold14.py`/`step_06b_gold14.py` integrate perim over `[π/2+1e-9, π/2+π/7−1e-9]`, deleting 1e-9 slivers at the **vertex** endpoints. Near a vertex the p-ball integrand has a boundary layer (|cos(θ−π/2)|ᵖ with θ−π/2∈[0,1e-9] amplifies each p-derivative by ≈ln(1e-9)=−20.7), so the deleted sliver contributes err(p) with err‴(1) ≈ +1.05e-5. Being **systematic** (the same smooth function at every p-sample), both Fourier radii and the FD stencil extract the *same* err‴(1) — hence self-consistent to 7e-18 yet all wrong by 1.05e-5. Verified by direct reproduction: gold-chamber(slivers) → −5.5080732424; full/centered chambers → the truth.
- **Control (12-gon, known P‴(1)=−5.2032404295176662):** the corrected pipeline (chamber *split* at the kink θ=0 — the 12-gon has a normal at π/2 so |sin θ|ᵖ vanishes at chamber center; no slivers) gives P(1)=12√(26−15√3) to 5e-70, Fourier P‴(1) to ~20 digits, FD to 2e-11. (The 14-gon has no normal at π/2, so its chamber center is smooth.)

---

## Task B — δ(N₃): **NONZERO** (like the octagon and 12-gon; not zero like the decagon)

I ported the 12-gon engine to **F = Q(2cos(2π/7))** exact arithmetic (`engine14.py`: FT field arithmetic, universal Jloc/Jtab tables, pair M1–M4/MixedW + TRIPLE recursions; validated vs quadrature to ~1e-59 on 13 integral types). Assembly: E3n14 → 811 reduced monomials → S (132 terms); **84·L·S = P‴(1) to 4.4e-50** (end-to-end validation).

σ(S) (weight-3 non-local symbol) has 8 terms, all in (span{α,u})³. Atom symbols: Li₃ = z⊗z⊗(1−z); A(x)=Li₃(x)−Li₃(1−x)+ζ(3)+prod (derived exactly, validated 1e-61); T3p via PGEN1R's repair (validated 1e-61); T3a via the S3 recursion on verified G-words; Lam via DOD3's S_Li21 (−τ into PGEN1R's convention).

- **Alphabet rule at the 14-gon (s=5):** σ(S) reduces to the 24 cross-ratio Li₃ symbols mod shuffle-decomposables with **EXACT residual 0** in all three F-components (Li₃-span mod Decomp is 2-dimensional). (EXACT)
- **δ(N₃) is well-defined:** it vanishes on all **22** independent Li₃ relations; identical across 4+ particular solutions, a minimal 2-element basis, and both Lam-sign conventions. (EXACT)

$$\delta(N_3)\;=\;c_\alpha\,(\alpha\wedge u)\otimes\alpha\;+\;c_u\,(\alpha\wedge u)\otimes u,\qquad \alpha=2\cos\tfrac{2\pi}{7},\ u=\alpha^2+\alpha-1=2\cos\tfrac{\pi}{7}$$

- **c_α = (601 − 182α − 243α²)/32 ≈ −0.118909811681569**
- **c_u = (613 − 118α − 292α²)/32 ≈ +0.369019755994578**

**δ(N₃) ≠ 0** — the 14-gon's far kink lines carry **genuine depth-2 (Li₂,₁) trilogarithms**, like the octagon and 12-gon, unlike the decagon. Derivation status: **EXACT / PROVED at symbol level** (exact engine, exact reduction, relation-invariant δ; atom identities proved + numerically cross-checked to 1e-61).

- **C₇ check (paper's local law):** 84·L·ζ(3)_local = **C₇ = −0.6838425376363484** to **3.6e-69** (C_m = 48m·sin⁴(π/2m)cos(π/2m)(3sin²(π/2m)−1), m=7). The ζ(3) coefficient in S is *entirely* local (far families contribute exactly 0). (EXACT)

---

## Task C — closed form of P‴(1): **OPEN**

Not attempted (δ(N₃) was the goal and is complete). The symbol-level reduction σ(S)=Σ c_zLi₃(z) mod products is in hand (`delta14.pkl`, `q3setup.pkl`); the full DOD5-style closure (products + invisible constants) remains open.

---

## Evidence table (all printed by code)

| check | result | label |
|---|---|---|
| Richardson FD P‴(1), centered chamber | −5.508083698053478284 | NUMERICAL 24d |
| 84·L·∫E₃n (exact chain) | −5.508083698053478284 | NUMERICAL 24d |
| gold chamber (slivers) → −5.5080732424 | bug located | NUMERICAL |
| 12-gon control (P(1)/Fourier/FD) | 5e-70 / ~20d / 2e-11 | NUMERICAL |
| engine14 vs quadrature (13 types) | ≤1.2e-58 | NUMERICAL |
| 84·L·S vs P‴(1) | 4.4e-50 | NUMERICAL |
| σ(S)→24-CR Li₃ mod Decomp | residual exactly 0 | EXACT |
| δ vanishes on 22 Li₃ relations | well-defined | EXACT |
| δ(N₃) = c_α(α∧u)⊗α + c_u(α∧u)⊗u | ≠ 0 (coeffs above) | EXACT |
| δ robust (minimal basis / both conventions) | match / both nonzero | EXACT |
| 84·L·ζ(3)_local = C₇ | 3.6e-69 | EXACT |

## Step files that produced these
`step_A1–A6` (P‴ resolution + 12-gon control), `step_B1–B17` (mons14, engine14.py, S assembly, σ, reduction, δ), `step_C1–C10` (Ax/T3p validations, C₇, exact residual, δ finalization + robustness), `step_V1–V2` (verify runs). Key artifacts: `engine14.py`, `sym14.py`, `mons14.pkl`, `S14.pkl`, `sigma14.pkl`, `delta_FINAL.pkl`, **`verify_pg6.py`** (PASS), `result.md`, `summary.md`, `progress.log`.

**Bottom line.** P‴(1) of the 14-gon is **−5.5080836980534783** (the exact-integrand chain; the polar gold was corrupted by a 1e-9 vertex-sliver deletion that systematically shifted err‴(1) by +1.05e-5). Its (2,1) invariant **δ(N₃) is NONZERO** — c_α(α∧u)⊗α + c_u(α∧u)⊗u with c_α=(601−182α−243α²)/32, c_u=(613−118α−292α²)/32 — so the 14-gon, like the octagon and 12-gon, has genuine depth-2 trilogarithm content; and its local ζ(3) term equals the paper's C₇ exactly.