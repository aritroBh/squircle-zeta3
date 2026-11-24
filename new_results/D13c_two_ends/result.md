# D13c — Two ends, one set of numbers: ±iπ/2 at p = 1 and at p = ∞

All numbers below were printed by the code listed at the end (run in this directory with the project venv). Labels: **PROVED / NUMERICAL / CONJECTURE**.

---

## 0. Groundwork (checks of the identity and index shift)

- The 450-digit `f_0…f_180` (`../D11b2_stokes_jump/step2_num450_N180.json`, key `f450`) match the exact D11 ζ-forms: residuals **0** at 80 dps for f₀…f₅ (`step_02.py`).
- The identity **`M_k = (√2/2)[ε^k]F(ε-1)`** (equivalently `F(η) = √2/(1+η) + √2 Σ_{k≥0} M_k(1+η)^k`) checked against direct quadrature of `F(η) = ∫₀¹(t(1-t))^{η/2}√cosh(ηl)dt` (`step_02.py`):

| η | residual |
|---|---|
| −0.2 | 6.03e‑57 |
| −0.3 | 9.58e‑54 |
| −0.5 | 1.1e‑37 |
| −0.7 | 1.69e‑22 |
| −0.9 | 4.67e‑7 (= truncation of the k ≤ 400 table near the convergence radius) |

and F(−0.2) = 1.27233782076052… reproduces D11's reference value. **NUMERICAL ✓** (identity and index shift confirmed).

---

## 1. The dictionary (task 1)

**Setup.** `M_k = (√2/2)[ε^k]F(ε-1)`, and the singularity of `F(ε-1)` at ε = 1 (η = 0) is the Gevrey-1 object of D11. A trans-monomial `T(η) = S η^b e^{-A/η}` in F, with u = 1-ε = -η, contributes to `[ε^k]` through the essential-singularity saddle of `u^b e^{A/u}(1-u)^{-N}` (N = k+1) at `u = ±√(A/N)`; the decaying branch (`2√(AN) = -√(πN)(1±i)` for A = ±iπ/2) gives, by a standard Gaussian transfer,

`[ε^k] u^b e^{A/u} ≈ (1/(2i√π)) A^{(2b+1)/4} e^{A/2} N^{-(2b+3)/4} e^{2√(AN)} (1+O(N^{-1/2}))`.

**The N-power `-(2b+3)/4 = -1` (D13b) forces `b = 1/2`.** Then requiring the per-branch constant to equal D13b's `κ_± = -√π/2` (both branches, summing to the VERIFIED law `-(√π/N)e^{-θ}cos θ`, θ = √(πN)) fixes the trans-monomial uniquely:

**D13b forces:** `b = 1/2`, `S = 2i√π` for A = +iπ/2 (and S = 2i√π for the conjugate; |S| = 2√π per branch), i.e.

`F(η) ⊃ 2i√π η^{1/2} e^{-iπ/(2η)} + 2i√π η^{1/2} e^{+iπ/(2η)}`,

which assembles to `κ_s = (is/√2) z_s^{1/2} e^{z_s/2} = -√π/2` — matching D13b's own saddle formula **symbolically exactly**. Amplitude √π/N ✓, phase π ✓.

**κ1 check** (`step_11.py`): the brief's `κ1 = [(24+12π+π²)+i(24-12π+π²)]/(48√π)` equals `conj(β_+)`, `β_s = 1/(2√z_s) + √z_s/2 - z_s^{3/2}/12`, to **1.1e-81**; `c_r = 0.841215087768`, `c_i = 0.0450118376852` ✓. (The `c3 = -1/2` input lives in the integrand-level germ B(z); the Borel germ B_F has a different O(ζ) term on the first sheet — §2 — which changes nothing at leading order; this plausibly feeds D13b's residual 2% tension in c_r.)

**Recomputing D11's cut integral with the correct Borel–Laplace prefactor** `F(η) = (1/η)∫₀^∞ e^{-t/η}B(t)dt` (this 1/η is mandatory when `B(t) = Σ f_n t^n/n!`): with `Disc B` from `c_+ = -√(2π)e^{3iπ/4}`,

`J(η) = (1/η)e^{-A/η}·2ic_+A^{-1/2}Γ(3/2)η^{3/2} = 2√π η^{1/2} e^{-A/η}` (constant printed: 3.5449077 = 2√π exactly).

At η = iε: `J = 2√π e^{3iπ/4} ε^{1/2} e^{-π/(2ε)}` — **same constant 2√π, same phase e^{3iπ/4}, same rate as D11, but power 1/2, not 3/2** (ratio printed: D11 = η × correct jump, the `e^{5iπ/4}` factor of the η^{1/2}-branch).

**Verdict (PROVED at the level of the saddle/dictionary mechanism, NUMERICAL confirmed):** D11's `2√π e^{3iπ/4} ε^{3/2} e^{-π/(2ε)}` is **right up to normalisation** — specifically it is missing the `1/η` of the Borel–Laplace transform `F = (1/η)∫e^{-t/η}B dt`; the correct jump is `2√π e^{3iπ/4} ε^{1/2} e^{-π/(2ε)}`. The exponent is **b = 1/2** (planner's heuristic right), the Stokes constant magnitude **2√π** (D11 right). Both ends agree: `|S| = 2√π`, `b = 1/2`, `A = ±iπ/2`.

---

## 2. Independent check: Borel transform from f_n (task 2)

From the 450-digit data, `b_n = f_n/n!`, singular fit at t = iπ/2 (`steps_04–07.py`):

| quantity | method | value | target | residual |
|---|---|---|---|---|
| branch exponent γ+1 | period-4 ratio e_n, per-class 3-param fits [100,172] | 1.4999966 ± 3.4e-6 | 3/2 | 3.4e-6 |
| amplitude `b_n(π/2)^n n^{3/2}/(√2cos(π(3-2n)/4))` | quadrisected deg-3 fits [120,180] | 1.00000000 ± 1e-7 | 1 | 1e-7 (rms 5e-13) |
| `c_+` in `B_F ~ c_+(1-t/A)^{1/2}` | re-derived | `-√(2π)e^{3iπ/4} = √π(1-i)` | D11's value | exact |
| local const. `B_F(A+ζ)/ζ^{1/2}` | from c_+(-1/A)^{1/2} | 2.0 (−2e-121 i) | 2 (= √2·D13b's s√2 after the √2/2 map) | exact |
| next term d1 in `2ζ^{1/2}(1+d1 ζ+…)` | c1(ρ) pattern (−0.803096 / +1.553038) | **d1 = (1-iπ)/4** | — | matches to 4 digits |

So the Borel transform has a **square-root branch (γ = 1/2)** at ±iπ/2 with local constant `B_F(A+ζ) ~ 2ζ^{1/2}` — NUMERICAL to 7e-6 (γ) and 1e-7 (amplitude). Via `J = 2iγ_AΓ(3/2)η^{1/2}e^{-A/η}` this gives exactly the D13b-forced `b = 1/2`, `|S| = 2√π`. **The two ends check each other: consistent.**

---

## 3. Next singularities (task 3)

**Prediction (PROVED as a mechanism, NUMERICAL confirmed):** writing `F(η) = 2∫₀^∞ e^{-l(1+η/2)}(1+e^{-l})^{-2-η}(cosh ηl)^{1/2}dl` and expanding `(1+e^{-l})^{-2-η} = Σ_j c_j(η)e^{-jl}` before substituting x = ηl, term j is a Laplace integral of `h(x) = (1+e^{-2x})^{1/2}` (branch points `x = iπ(2m+1)/2`, local constant √2) against `e^{-x(1+j)/η}`; hence **F has Borel singularities at t = ±iπq/2 for every q = 1,2,3,…** (`q = (2m+1)(1+j)`), each of square-root type, giving trans-monomials `η^{1/2}e^{-A_q/η}` and, through the same saddle, `e^{-√(qπN)}` pieces in M_k. So:

- `e^{-√(2πN)}` ⟸ **A = ±iπ** (q = 2: the m = 0, j = 1 tail term), predicted `R_k ≈ (√π/N)e^{-√(2πN)}cos(√(2πN)-π/4)`;
- `e^{-√(3πN)}` ⟸ **A = ±3iπ/2** (q = 3: both (m=1,j=0) and (m=0,j=2), predicted S₃ = 2i√π(1±1/√3), sheet sign unfixed).

**Checks** (`steps_08–10.py`; M_k^main evaluated by rotated-contour quadrature, cached `mkmain_cache.json`):

- The exact 1D integral `M_k^main = (2/π)Im∫₀^∞ e^{iπs/2}(s/(1+s))^N H(s)ds/s` matches the 52-digit table to O(e^{-√(2πN)}): e.g. k=100: diff 1.7e-13; k=400: diff 4.8e-25 (diff/scale ∈ [0.001,0.016]) — **independent verification of D13b's exact integral**.
- Rate test on R_k = M_k − M_k^main (k = 100…400): **a = √(2π) wins decisively** (fit rms 1.4e-19 vs 7.2e-18 (√(3π)), 1.3e-17 (2√π), 1.0e-17 (√π); wrong rates force absurd amplitudes 2.5e3, 9.3e5).
- With the N^{-1/2} correction included, the fit `R_k = e^{-√(2πN)}N^{-1}[(C cos ω2 + D sin ω2)+…]`, ω2 = √(2πN), over windows [kmin,400]: C = 1.264→1.252, D = 1.140→1.226, amp = 1.702→1.752, phase = 0.734→0.775 — trending to the prediction **C = D = √(π/2) = 1.2533, amp = √π = 1.7725, phase = π/4 = 0.7854**. Frequency detune scan: minimum at δ ≈ +0.005 (consistent with the N^{-1/2} phase drift). **NUMERICAL ✓ — the e^{-√(2πN)} piece is confirmed, including its amplitude and phase, from the F-side.**
- q = 3 (e^{-√(3πN)}): after subtracting the fitted q = 2 piece, the residual (≲1e-24 at k ≥ 200) is comparable to the q = 2 model error; **not resolved** at k ≤ 400. Prediction recorded above; status **CONJECTURE (untested)**.

---

## 4. Rigour for D13b (task 4)

Working from the exact 1D integral (independently re-verified above):

- **Contour rotation** to the ray arg s = π/4: justified by `|H(s)| ≤ H(Re s)` (H(s) = ½∫₀¹t^{s/2-1}(1-t)^{1/2}dt) and the sectorial decay of `e^{iπs/2}`. **PROVED.**
- On the ray, `Re Φ(x) = -πx/(2√2) + N log x - (N/2)log(1+√2x+x²)` is **strictly concave for all x > 0**: `Re Φ''(x) = -N/x² + Nx(x+√2)/(1+√2x+x²)² < 0` ⟺ `(1+√2x+x²)² - x³(x+√2) = √2x³ + 4x² + 2√2x + 1 > 0` — a polynomial identity. **PROVED** (also scanned: 0 sign violations on 4000 points).
- Exact saddle `s* = -½ + √(¼ + 2iN/π)`; sympy: `Φ(s*) = -√(πN)(1-i) - iπ/4 + O(N^{-1/2})`, `Φ''(s*) = √2√i·iπ^{3/2}/(2√N)`. **PROVED (symbolic).**
- All Laplace constants are explicit and N-stable (`step_14.py`): `(ReΦ(x0)+θ)√N → 0.44`, `|Φ''(x0)|√N → 3.7`, `|ReΦ''(x0)|√N → 2.78`, curvature `c(N) ≈ 0.414/√N`, `|ψ(x0)|N^{5/4} ≈ 2.1`, `|dlogψ| ≤ 0.31`. The total error ratio `M_k^main/(-(√π/N)e^{-θ}cosθ) - 1 = -0.122, -0.075, -0.060, -0.046` at k = 50,100,200,400 — O(1/√N) as required. **NUMERICAL.**

**Status: NOT yet PROVED — one step is missing.** What is proved: exact representation, contour rotation, concavity (so a unique saddle and global two-sided bounds `ReΦ(x) ≤ ReΦ(x0) - c(N)(x-x0)²` on `[x0/2, 3x0/2]` plus linear tails), and the symbolic saddle value whose `(2/π)Im` equals `-(√π/N)e^{-θ}cosθ(1+O(N^{-1/2}))`. **Missing:** the final assembled stationary-phase remainder inequality — choosing the window W = σw_N, bounding `|R₃| ≤ M₃|δ|³/6`, `|ψ-ψ(x0)|`, and the tail integrals with the explicit constants above, to produce a printed constant C and threshold N0 in `|M_k^main + (√π/N)e^{-θ}cosθ| ≤ C N^{-3/2}e^{-θ}` — the computation of that constant was interrupted (all ingredients are printed and bounded, but C and N0 were never printed, so I do not claim them). The `m ≥ 1, n ≥ 1` pieces are bounded **NUMERICALLY** by the e^{-√(2πN)} residual fits of §3 (their rigorous bound requires the full M(q) representation and is not assembled here).

---

## What this means for the squircle

The same two numbers ±iπ/2 control both ends. At p = 1 they are the zeros of cosh z, which put square-root branch points at ±iπ/2 in the Borel plane of F, make the rhombus Taylor series diverge factorially (f_n ~ √2 n!(2/π)^n n^{-3/2}), and set the Stokes jump `2√π η^{1/2} e^{∓iπ/(2η)}`. At p = ∞ the very same pair — now as the branch points of the integrand-level B(z) — become the essential singularity of F(ε-1) at ε = 1, whose saddles turn the exponential e^{∓iπ/(2η)} into e^{-√(πN)}cos(√(πN)): the squircle's large-k coefficients are the Borel plane of its rhombus end, read through one saddle. One local constant (the square-root branch of magnitude 2√π at ±iπ/2) fixes both the Stokes jump at p = 1 and the decay law at p = ∞; D11's jump was right except for one missing power of η — the Laplace 1/η — which this cross-check exposed.

---

## Evidence table

| check | result | residual / error | digits |
|---|---|---|---|
| f_n vs exact ζ forms (n ≤ 5) | match | 0 @80dps | 80 |
| M_k = (√2/2)[ε^k]F(ε-1) vs quadrature, η ∈ (-0.7,-0.2) | ✓ | 1.7e-22 … 6e-57 | 22–57 |
| κ1 (brief) = conj(β_+) | ✓ | 1.1e-81 | 81 |
| Borel exponent γ | 0.4999966 ± 3.4e-6 | 7e-6 | ~6 |
| Borel amplitude | 1.00000000 ± 1e-7 | 1e-7 | 7 |
| γ_A = B_F(A+ζ)/ζ^{1/2} | 2.0 (−2e-121 i) | 1e-121 | 121 |
| correct jump constant (with 1/η) | 2√π | printed 3.5449077 = 2√π | exact |
| M_k^main vs M_k | diff = O(e^{-√(2πN)}) | 1.7e-13 (k=100) … 4.8e-25 (k=400) | 13–25 |
| q=2 rate test | a = √(2π) wins | rms 1.4e-19 vs ≥7.2e-18 others | — |
| q=2 amplitude/phase fit | → √π = 1.7725, π/4 = 0.7854 | 1.752 / 0.775 at [300,400] | ~2 |
| ReΦ concavity | `√2x³+4x²+2√2x+1 > 0` | PROVED + 0/4000 violations | — |
| leading-law ratio M_k^main/law − 1 | −0.046 at k=400, O(1/√N) | NUMERICAL | — |

## Files (in this directory)

- `step_02.py` — data loading; f_n vs ζ forms; identity M_k = (√2/2)[ε^k]F(ε-1) vs quadrature.
- `step_03.py` — bare-model cautionary check (growing saddle; motivates branch-cut transfer).
- `step_04.py`, `step_05.py`, `step_06.py`, `step_07.py` — Borel transform fits: exponent γ (1.4999966 ± 3.4e-6), amplitude (1 ± 1e-7), c_+, γ_A = 2, d1 = (1-iπ)/4.
- `step_08.py` — exact 1D integral M_k^main by rotated quadrature; matches M_k to O(e^{-√(2πN)}); cache `mkmain_cache.json`.
- `step_09.py`, `step_10.py` — residual fits: rate √(2π) decisive; C,D,amp,phase → √(π/2), √π, π/4; frequency detune.
- `step_11.py` — q=3 attempt (unresolved); κ1 = conj(β_+) to 1e-81; D11 jump recomputation with the 1/η (constant 2√π, power 1/2).
- `step_12.py`, `step_13.py`, `step_14.py`, `step_15.py` — rigour: exact saddle and phase series (sympy), concavity proof + scan, Laplace constants table, leading-law ratios.
- `step_16.py` — saddle-value assembly and κ-fit of M_k^main (interrupted; see §4 missing step).