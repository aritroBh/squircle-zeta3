md = r'''# ZLAD1R — referee report on ZLAD1 (the ζ-ladder at every angle; universal Γ-kernel)

Work dir: `research/08_theory/ZLAD1R_referee_ladder/`. All computations are **independent re-derivations**
(no ZLAD1 code used): my own exact table `muR_exact_N40_K8.json` is built from `thm:flat`'s Beta-integral
formula via a different route than ZLAD1 (binomial expansion of `f^{1-2n}` with Stirling numbers +
direct bivariate Γ-series exponentiation, instead of ZLAD1's Touchard polynomials + directional Bell
recursion), and it also differs from G1's `flat_exact.py` (numeric n substituted first; no symbolic-n
sympy series). Cross-checks: G1 (k≤3, full range n≤300), G2 (k=4, n≤300), ZLAD1 (k≤8, n≤150).

## 1. The objects

`M(u;α) = Σ_n μ(u;n) c^{2n}`, `μ(u;n) = -G_n(u)/(n(2n-1))`,
`G_n(u) = B(n,n)^{-1} ∫_0^1 y^{n-1}(1-y)^{n-1} E_n dy`, `E_n = e^{uA} f^{1-2n}`,
`A = -(n-1)a + nb`, `f = 1-y+y e^{uℓ}`, `ℓ = b-a`, `a = log y`, `b = log(1-y)`.
Integration: `y^j a^i b^l ↦ (n)_j/(2n)_j · D_{i,l}(n+j,n)`, `D_{i,l}` = Bell polynomials in the
cumulants `κ_{r,s} = ∂_p^r∂_q^s log B`. Zetas enter only through the mixed cumulants `κ_{r,s}, r,s ≥ 1`,
whose zeta part is `C_{r+s} ζ(r+s)`, `C_m = (-1)^{m+1}(m-1)!` (all `rs=0` cumulants: zeta and γ cancel).
All confirmed independently (step_22: γ cancels exactly in every `μ_{k,n}`, k≤8, n≤40 — 0 leftovers).

## 2. Task 1 — Lemma 1 and Lemma 2, line by line

**Lemma 1 (operator formula).** `[ζ^E]G_n(u) = [𝒫^E(∂ξ,∂η)Ψ_n](u;0,0)`,
`𝒫^E = ∏_w ((C_w/w!)□_w)^{e_w}/e_w!`, `□_w = 𝒟^w - ∂ξ^w - ∂η^w`, `𝒟 = ∂ξ+∂η`.
- `e^Δ = e^{Δ^H}e^{Δ^Z}` with `Δ^Z = Σ_w ζ(w)(C_w/w!)□_w(x,y)`: **CORRECT.** `Δ^Z` is ζ-linear, so
  `[ζ^E]e^Δ = e^{Δ^H}·𝒫^E(x,y)` as formal series (polynomial in ζ's, finite at each degree).
- `[x^iy^l]`-extraction = applying `𝒫^E(∂a,∂b)` before the 𝕳-map: **CORRECT** (Leibniz: both sides are
  `Σ_{r,s} p_{rs} (i)_r(l)_s D^H_{i-r,l-s}`; the `(p,q)=(n+j,n)`-dependence is untouched since 𝒫 acts
  only on `a,b`). Verified numerically with MY code at 17 points (n,k,E) incl. ζ(7), ζ(8), ζ(2)⁴,
  ζ(3)ζ(5), ζ(2)²ζ(3), ζ(2)ζ(5), ζ(4)², ζ(2)ζ(6): **ALL EXACT** (step_28).
- Reduction modulo `𝒟-u` (used next): valid because `(𝒟-u)Ψ_n = 0` **exactly** (see Lemma 2);
  operators are commutative polynomials, so `Q·(𝒟-u)` kills Ψ_n for any Q. **CORRECT.**
- Step back to `μ_{k,n}`: coefficient-wise in `u` with `c_e^{(E)}(u) ∈ ℚ[u]`: **CORRECT.**

**Lemma 2.** `(∂a+∂b)E_n = uE_n` exactly: `𝒟A = -(n-1)+n = 1`, `𝒟f = 0` since `f` depends on `(a,b)`
only through `ℓ = b-a`. **CORRECT** (verified with MY code, `n∈{2,3,5}`, `k≤8`: step_26).
The stated isomorphism `ℚ(u)[∂ξ,∂η]/(𝒟-u) ≅ ℚ(u)[s]` in the summary result.md is **MINOR FIX**
(the quotient is `ℚ(u)[ξ]`, one variable; the *symmetric* operators land in `ℚ(u)[s]`, `s = ξ(u-ξ)` —
and `□_w` IS symmetric, so nothing is affected; the long writeup `step_44.py` states it correctly).
There `□_w ↦ β_w(u,s) = u^w - p_w(u,s)`, `p_w` Lucas (`p_0=2, p_1=u, p_w = u p_{w-1} - s p_{w-2}`):
**CORRECT** (β_2 = 2s, β_3 = 3us, β_4 = 4u²s - 2s²; reproduced independently).

## 3. Task 2 — Theorem ZLAD1′ (universal kernel)

`Σ_w (-1)^{w+1} ζ(w) β_w/w = log[Γ(1+ξ)Γ(1+η)/Γ(1+u)]`, `ξ+η = u`, `ξη = s`:
- γ-cancellation: `logΓ(1+t) = -γt + Σ_{m≥2}(-1)^m ζ(m)t^m/m`; the `-γ` parts sum to `-γ(ξ+η-u) = 0`.
  **CORRECT.** (Also: γ cancels identically inside every mixed Beta kernel — step_22.)
- The identity itself: **CORRECT** — verified symbolically EXACT to `u^13` (lhs - rhs ≡ 0, step_30)
  and numerically at 5 `(u,s)` points inside the polydisc, worst residual 8.4e-39 (step_31).
- Resummation `G_n(u) = Σ_{e≥0} T_e Φ_e`, `Φ_e = (-∂_s)^e Ω(u,0)`, `Ω = Γ(1+ξ)Γ(1+η)/Γ(1+u)`:
  follows from ZLAD1 since `Σ_E ζ^E 𝒫^E ↦ exp(Σ_w ζ(w)(C_w/w!)β_w) = Ω` (`C_w/w! = (-1)^{w+1}/w`).
  **CORRECT** — verified EXACTLY (not just numerically) for all `n=2..40`, `k≤8`, every zeta monomial
  (2609 values, step_30) with MY table. `Φ_0 = 1`, `Φ_1 = (ψ(1+u)-ψ(1))/u = Σ_{w≥2}(-1)^w ζ(w) u^{w-2}`,
  `Φ_2 = Φ_1² + Σ_{w≥4}(-1)^w (w-3) ζ(w) u^{w-4}` reproduced.

## 4. Task 3 — independent tests

- **Table** `μ_{k,n}`, `k≤8`, `n≤40`, from the Beta-integral formula (own code, step_22):
  == ZLAD1 `mu_exact_N150_K8.json` exactly (all 2676 monomial coefficients, step_23);
  == G1 (k≤3) and G2 (k=4) exactly; ZLAD1's K6 table == G1 (n≤300) and G2 (n≤300) exactly (step_32).
- **Quadrature end-to-end** (Cauchy extraction of `[u^k]` of the actual integral, 60 digits):
  (k,n) = (8,3), (7,5), (5,7): residuals 3.8e-60, 5.6e-61, 4.7e-62 (step_24).
- **Five rungs** (EXACT, my table, all n≤40; step_25):
  `[ζ(3)]μ_8 = -[ζ(2)]μ_7`; `[ζ(7)]μ_8 = -[ζ(2)]μ_3 -4[ζ(2)²]μ_5 -6[ζ(2)³]μ_7` (weight 7);
  `[ζ(8)]μ_8 = [ζ(2)]μ_2 +5[ζ(2)²]μ_4 +12[ζ(2)³]μ_6 +6[ζ(2)⁴]μ_8` (weight 8);
  `[ζ(3)ζ(5)]μ_8 = 2[ζ(2)²]μ_4 +6[ζ(2)³]μ_6` (mixed, w 8); `[ζ(2)²ζ(3)]μ_7 = -3[ζ(2)³]μ_6` (mixed, w 7).
  Plus the remaining 12 mixed/single rungs at k=7,8 (step_26): all PASS.
- **Closed forms** (EXACT, n≤40; step_25): `[ζ(5)]μ_{5,n} = -3/((2n-1)(2n+1)(2n+3))` and
  `[ζ(2)³]μ_{6,n} = n(5n-1)/(3(2n-1)(2n+1)(2n+3)(2n+5))` (the two requested), plus
  `[ζ(2)ζ(3)]μ_5`, `[ζ(3)²]μ_6`, `[ζ(6)]μ_6`, `[ζ(2)ζ(4)]μ_6`, `[ζ(2)⁴]μ_8`, `[ζ(8)]μ_8`, `[ζ(4)²]μ_8`:
  all 9 PASS.
- Function level (c²=1/2): rung 1 holds for k=3..8; `[ζ(2)]M_2(π/2) = Φ(1/√2)` to 1.3e-16 (step_33).

## 5. Verdict

| Item | Verdict |
|---|---|
| Lemma 1 (operator formula) | **PROVED** (all k, n); steps verified; 17-point independent exact check |
| Lemma 2 (`(∂a+∂b)E_n = uE_n`; reduction) | **PROVED**; summary's `≅ ℚ(u)[s]` = MINOR FIX (cosmetic) |
| Theorem ZLAD1 (master ladder) | **PROVED** (all k, n) |
| Theorem ZLAD1′ (universal Γ-kernel) | **PROVED**; γ-cancellation and logΓ identity correct; resummation verified EXACTLY n≤40, k≤8 |
| Rung list (6 single + 10 mixed) | **PROVED**; all verified EXACTLY at k≤8, n≤40 on an independent table |
| Closed forms `[ζ(5)]μ_5`, `[ζ(2)ζ(3)]μ_5`, `[ζ(3)²]μ_6` | **PROVED** (inputs: thm:flat's `1/(4n²-1)` + G2's proved `-n/((2n+3)(4n²-1))`); re-verified n≤40 |
| Free-tower closed forms `[ζ(2)³]μ_6`, `[ζ(2)⁴]μ_8` and descendants `[ζ(6)]μ_6`, `[ζ(2)ζ(4)]μ_6`, `[ζ(8)]μ_8`, `[ζ(4)²]μ_8` | **EXACT (not proved for all n)**: verified n≤40 (mine), n≤300/150 (theirs) — labeling in result.md is honest |
| All-k ladder `[ζ(3)]M_{k+1} = -[ζ(2)]M_k` | **PROVED** (rung 1 of ZLAD1) |

**Paper statement (proposal).** *Theorem (ζ-ladder). Let M(u;α) = Σ_n μ(u;n) c^{2n} be the flat-angle
expansion of Theorem thm:flat. Then every zeta-monomial coefficient of μ(u;n) is a universal
ℚ[u]-combination of the pure ζ(2)-power coefficients: [ζ^E]μ(u;n) = Σ_{e≥1} (-1)^e e! c_e^{(E)}(u)
[ζ(2)^e]μ(u;n), where ∏_w ((C_w/w!)β_w)^{e_w}/e_w! = Σ_e c_e^{(E)}(u) s^e with
β_w = u^w - p_w(u,s) (p_w Lucas) and C_w = (-1)^{w+1}(w-1)!. Equivalently, the zeta values enter
M(u;α) only through the single universal kernel Ω(u,s;ζ) = Γ(1+ξ)Γ(1+η)/Γ(1+u) (ξ+η = u, ξη = s):
M(u;α) = Σ_{e≥0} 𝒯_e(u;α) (-∂_s)^e Ω(u,0;ζ), 𝒯_e := [ζ(2)^e]M(u;α). In particular
[ζ(3)]M_{k+1}(α) = -[ζ(2)]M_k(α) for every k ≥ 0 at every angle.*

## 6. Files

`step_22.py` (independent table), `step_23.py` (cross-diffs), `step_24.py` (quadrature),
`step_25.py`/`step_26.py` (rungs + closed forms), `step_27.py`/`step_28.py` (Lemma 1),
`step_26.py` (Lemma 2), `step_29.py`/`step_30.py`/`step_31.py` (ZLAD1′),
`step_32.py` (full-range diffs), `step_33.py` (π/2), `verify_zlad1r.py` (< 5 min).
'''
open('result.md', 'w').write(md)
open('progress.log','a').write('step_35: result.md written (early); all main checks done\n')
print('result.md written', len(md))
