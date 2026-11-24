# D15c — The edge law for ζ(3) near p = 1, and C_m in closed form

**Directory:** `research/08_theory/D15c_edge_zeta3/`. Predecessor `../D15b_polygon_law/` (result.md, poly_coeffs.py, step_09.py, step_10.py, z3coeff_m2_8.json, vals_m2_8_200.json) read first; `../R2_rhombus_ab/` read for the proof technique and the rhombus target.

Labels: **PROVED** = exact symbolic algebra or proved consequence; **NUMERICAL** = high-precision agreement with printed residual; **CONJECTURE** = pattern without proof.

---

## 1. Task (1) — Symbolic G₂₁, G₂₁₂: PROVED + validated to 1e-61

On the edge arc only s₀=sin θ and s₁=s_{m−1}=sin(A−θ) vanish. Writing L₀=ℓ₀, L₁=ℓ₁ as symbols and separating the **reduced** non‑vanishing sums τₖ=Σ_{i∈N}sᵢℓᵢᵏ, γₖ=Σ_{i∈N}cᵢ(ℓᵢᵏ+kℓᵢᵏ⁻¹), I built D₃=D₀(ρ₃+ρ₂E₁+ρ₁E₂+E₃) exactly in sympy (`sym_G.py`), expanded, and extracted the coefficients of L₀²L₁ and L₀L₁². Because S₀=T₀=Σsᵢ and B=1+Q₀² are log‑free, D₃ is a **polynomial of total degree ≤3** in (L₀,L₁), so step_09's degree‑3 perturbative fit returns the *true* coefficients.

**Structural theorem (PROVED, `step_03_check.py`):** G₂₁ and G₁₂ depend **only on s₀,s₁,c₀,c₁,S₀,Sc** — they have *no* τₖ, γₖ, or log S₀ dependence. So the "non‑vanishing logs" never enter the cross coefficients; **the Clausen‑type terms D15b feared are absent from G₂₁,G₁₂ identically** (not just cancelling on integration). Closed form (validated below):

G₂₁ = −√(S₀²+Sc²)·N₂₁ /(2 S₀⁵(S₀²+Sc²)^{5/2}), with N₂₁ the explicit degree‑9 polynomial printed by `step_04_factor.py`; G₁₂ = G₂₁ with (s₀,c₀)↔(s₁,c₁) (the θ↔A−θ symmetry, verified `step_05_subst.py`: G₁₂(φ)−G₂₁(−φ)=0 exactly).

**Validation (NUMERICAL, `step_01_validate_G.py`, 60 dps):** symbolic G₂₁,G₁₂ vs D15b's perturbative extractor at 5 points (m=2,3,4,5): max residual **5.1e‑61** (G₂₁), **7.6e‑61** (G₁₂) — better than the required 50 digits.

## 2. Task (2) — C_m in closed form: PROVED

With φ=θ−A/2, h=A/2=π/(2m): S₀=cosφ/sin h, Sc=−sinφ/sin h, hence **S₀²+Sc²=1/sin²h (constant)** and √(1+Sc²/S₀²)=1/cosφ. Substituting, G₂₁+G₁₂ collapses (`step_07_integ2.py`) to **just 6 terms**, only sec²φ and sec⁴φ with coefficients in (sin h,cos h); the sf‑odd part vanishes identically (0 dropped terms). Integrating ∫sec²=tan, ∫sec⁴=tan+tan³/3 over φ∈[−h,h] (`step_08_integ3.py`):

$$\int_0^{A}(G_{21}+G_{12})\,d\theta \;=\; 2\,\sin^4 h\,\cos h\,\bigl(3\sin^2 h-1\bigr),\qquad h=\frac{\pi}{2m}.$$

$$\boxed{\;C_m \;=\; 24m\!\int_0^{A}\!(G_{21}+G_{12})\,d\theta \;=\; 48\,m\,\sin^4\!\tfrac{\pi}{2m}\,\cos\!\tfrac{\pi}{2m}\,\Bigl(3\sin^2\!\tfrac{\pi}{2m}-1\Bigr)\;}$$

Equivalent single‑radical form: C_m = 6m(1−cos(π/m))²(1−3cos(π/m))·√((1+cos(π/m))/2).

**Verification (NUMERICAL, `step_09_verify.py`, 130 dps):** formula vs `z3coeff_m2_8.json` for m=2..8, residuals **2.3e‑121 … 7.1e‑121**; vs known closed forms exactly: C₂=6√2 (resid 0), C₃=−9√3/8 (resid 0), C₆=(189√6−333√2)/8 (resid 7.7e‑130), C₈=−9g⁷+96g⁵−336g³+384g (resid 1.7e‑129). **Prediction + confirmation (`step_10_extract_new.py`, 100 dps):** running D15b's independent extractor for m=9,10,12 matches the formula to **7.1e‑102, 7.1e‑102, 3.6e‑102**.

The sign change of C_m is at 3sin²(π/2m)=1 ⟺ π/2m=π/6 ⟺ m=3 between m=2 (positive) and m=3 (negative): consistent with D15b's per‑corner sign change between α=π/2 and 2π/3; |C_m|→0 like −9π⁵/(16m⁴).

## 3. Task (3) — the non‑cross terms contribute no ζ(3)

Statement to prove: [ζ(3)]∫₀^A D₃dθ = 2ζ(3)∫₀^A(G₂₁+G₁₂)dθ (equivalently the reduced integrand R=D₃−cross is ζ(3)‑free).

* **Mechanism (PROVED).** The Beta moment B(α+u,β+v)/B has log = [logΓ(α+u)−logΓα]+[logΓ(β+v)−logΓβ]−[logΓ(α+β+u+v)−logΓ(α+β)]. The **pure‑direction** u³ coefficient is ψ″(α)−ψ″(α+β) = −2 at (1,1) — rational, **ζ(3)‑free by telescoping**; the **cross** u²v, uv² coefficients carry the universal +ζ(3). Thus [ζ(3)]∫g·ℓ₀²ℓ₁ = 2ζ(3)∫g and [ζ(3)]∫g·ℓ₀ℓ₁² = 2ζ(3)∫g for analytic g (linearity + uniform convergence), while ∫g·logᵃt·logᵇ(1−t) with (a,b)∉{(2,1),(1,2)} has no *universal* ζ(3). This is R2's ψ″(1)=−2ζ(3) localized to an edge.
* **The cancellation is real, not trivial.** Single‑corner log³ *does* carry ζ(3) in isolation: I printed ∫₀^{π/2}log³(sinθ)dθ = −(π/2)log³2−(π³/8)log2−**(3π/4)ζ(3)** (`step_16_settle2.py`, residual r₁/(πζ(3))=−0.75 exactly). c₃₀,c₀₃ are pure rational‑trig; c₂₀,c₀₂,c₁₁,c₀₀ carry the non‑vanishing logs (`step_11_allcoeff.py`). The non‑cross pieces therefore cancel only in sum.
* **Anchors (PROVED/NUMERICAL).** m=2: the independent full total is R2's P‴(1)=6√2ζ(3)−2√2π², and my `fullint−2ζ(3)·crossG = −√2π²/12` exactly (`step_14_rem.py`, and `step_17_pslq3.py` residual 0.0) → cross=total, PROVED. m=3: cross value −9√3/8 = independently planner‑checked total (D15b, 1e‑100). 
* **Independent cross‑check (NUMERICAL, `step_14_rem.py`):** 12m·fullint reproduces D15b's independent P‴(1) from `vals_m2_8_200.json` to **1.8e‑94 … 7.9e‑95** for m=2..8, and crossG=C_m/(24m) to **≤4.5e‑103**. Hence rem_m=fullint−2ζ(3)·crossG satisfies 12m·rem_m=P‴(1)−C_mζ(3); identifying rem_m's ζ(3)‑free closed form for general m (the remaining hard PSLQ) is the one piece I did not finish. **Status: cross‑term rule and mechanism PROVED; cross=total PROVED for m=2 and anchored for m=3; general cancellation NUMERICAL (≤1e‑94 via independent totals) with the telescoping structure identified.**

## 4. Task (4) — non‑regular test: rhombus (NUMERICAL, strong)

`step_18_rhombus.py`: weighted m=2 extractor for |x|^p+|λy|^p (semi‑axes a=1,b=1/λ). The cross‑term ζ(3) coefficient reproduces **R2's 24a²b²/r³** for λ = 1, 2, 3, 0.5, 1.7 with residuals **5.1e‑80 … 1.5e‑79**. This is the edge law on a non‑regular body: the rhombus's ζ(3) is fixed by its edge data (side length via r, corner angles π/2) exactly as R2 proved. The three‑different‑weights hexagon run was not completed before the turn limit (framework in place).

**Edge law (supported by all evidence):** the ζ(3) of P‴(1) is a sum over edges; each edge contributes only through the interference of its two endpoint corner logs (the ℓ₀²ℓ₁+ℓ₀ℓ₁² cross terms), with coefficient depending on that edge's data — for the regular 2m‑gon, on the corner angle α=π(m−1)/m alone, giving C_m = 48m sin⁴(π/2m)cos(π/2m)(3sin²(π/2m)−1).

## 5. Task (5) — "ζ(3) lives on the edges"

A polygon's boundary is flat almost everywhere: along each side nothing curves, and the curvature that builds the perimeter correction is concentrated at the corners. One might expect ζ(3) — the simplest odd zeta value — to be a "corner" number, stamped out at each vertex. It is not. A single corner, no matter its angle, produces only elementary constants: powers of π, logarithms, Clausen values. The number ζ(3) appears only where **two** corners talk to each other across the flat side that joins them. Each side of the polygon is bordered by a left corner and a right corner; the logarithm of the distance to each endpoint is a corner object, and it is precisely the *mixed* term — the square of one endpoint's log times the other endpoint's log — that carries ζ(3), with the universal coefficient 2ζ(3) inherited from the Gamma function's ψ″(1). Add up one such mixed contribution per side and you have the whole of ζ(3) in P‴(1). So ζ(3) is not a corner invariant and not a global one: it lives on the edges, in the interference between the two vertices that each edge connects.

---

## Evidence table

| Item | Check | Residual | Digits | Status |
|---|---|---|---|---|
| G₂₁,G₂₁ symbolic vs extractor, 5 pts | step_01 | ≤7.6e‑61 | 60 | PROVED + NUM |
| G₂₁,G₁₂ have no τ/γ/lS0 dep | step_03 | exact (0 monoms) | — | PROVED |
| G₁₂(φ)=G₂₁(−φ) | step_05 | exact 0 | — | PROVED |
| ∫(G₂₁+G₁₂)=2sin⁴h cosh(3sin²h−1) | step_07/08 | exact (6 terms) | — | PROVED |
| C_m formula vs json, m=2..8 | step_09_verify | ≤7.1e‑121 | 130 | PROVED + NUM |
| C₂=6√2, C₃=−9√3/8 | step_09_verify | 0.0 | 130 | PROVED |
| C₆,C₈ cyclotomic forms | step_09_verify | 7.7e‑130, 1.7e‑129 | 130 | NUM |
| m=9,10,12 prediction | step_10 extractor | ≤7.1e‑102 | 100 | NUM (confirmed) |
| 12m·fullint vs independent P‴ | step_14 | ≤7.9e‑95 | 100 | NUM |
| crossG = C_m/(24m) | step_14 | ≤4.5e‑103 | 100 | NUM |
| ∫₀^{π/2}log³sinθ ζ(3) coeff | step_16 | r₁/(πζ3)=−3/4 exact | 80 | NUM |
| m=2 full P‴=6√2ζ3−2√2π² | step_17 | 0.0 | 150 | PROVED (R2) |
| Rhombus cross = 24a²b²/r³, 5 λ | step_18 | ≤1.5e‑79 | 80 | NUM (edge law) |

## Files produced (in `D15c_edge_zeta3/`)
- `sym_G.py` — symbolic builder of D₃, G₂₁, G₁₂ (core).
- `step_01_validate_G.py` — validates symbolic G vs extractor (≤7.6e‑61).
- `step_02_simp.py`, `step_03_check.py`, `step_04_factor.py`, `step_05_subst.py` — simplification, τ/γ‑independence proof, factoring, symmetry.
- `step_06_integ.py`, `step_07_integ2.py`, `step_08_integ3.py` — the 6‑term collapse and exact integration → I(h); `monolist.pkl`, `I_of_h.txt`.
- `step_09_verify.py` — C_m formula vs json/known forms + predictions; `pred_m9/10/12.txt`.
- `step_10_extract_new.py` — independent extractor for m=9,10,12 (confirms formula).
- `step_11_allcoeff.py`, `allcoeffs.pkl` — all 10 log‑power coefficients of D₃.
- `step_12_explore.py`–`step_16_settle2.py` — decomposition, rem_m (`rem_m.json`), independent‑P‴ cross‑check, log³‑sine ζ(3) confirmation.
- `step_17_pslq3.py` — m=2 full P‴ exact.
- `step_18_rhombus.py` — rhombus → R2 edge law (≤1.5e‑79).

**Bottom line.** The ζ(3) part of P‴(1) for the regular 2m‑gon is **C_m·ζ(3) with C_m = 48m·sin⁴(π/2m)·cos(π/2m)·(3sin²(π/2m)−1)**, derived symbolically (G₂₁,G₁₂ exact, integrate to 2sin⁴h cos h(3sin²h−1), no Clausen terms — they vanish identically), verified to ~1e‑120 for m=2..8 and confirmed by the independent extractor for m=9,10,12. The mechanism (cross‑term Beta rule + ψ″ telescoping) is proved; cross=total is proved for m=2 (R2) and anchored for m=3; the general cancellation is established numerically to ~1e‑94 against independent totals. The rhombus reproduces R2's 24a²b²/r³ to 1e‑79, confirming the edge law on a non‑regular body. Open: the fully general symbolic proof of the non‑cross cancellation for all m (rem_m closed forms), and the three‑weights hexagon run.