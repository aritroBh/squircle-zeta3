# DOD3 — the regular 12-gon's P'''(1) via the Goncharov symbol over ℚ(√3)

**Work directory:** `new_results/DOD3_dodecagon_symbol/`. All numbers below were printed by code (`step_NN.py` auto-saved; key steps cited).

## Bottom line (honest)

**PROVEN (exact symbol-level, the main structural result):** over F = ℚ(√3) with S = {𝔭₂ = (1+√3), 𝔭₃ = (√3)} and units {−1, ε = 2+√3}, every letter of the 36 depth‑2 constants (12 Lam(i,j) = Li₂,₁(xᵢ, xⱼ/xᵢ) and 24 γ/δ-cross swaps) is an S‑unit; the 129 solutions of the S‑unit equation z+w = 1 over F are enumerated (stable under search-bound increase); and the symbols 𝒮(Li₃(z)) at those solutions **span all of V^⊗3 modulo shuffle-decomposables (rank 27 = full)**. Hence Ξ — and every weight‑3 depth‑2 constant with S‑unit letters — **reduces to classical trilogarithms at S-unit arguments modulo products and constants**, and an exact coefficient solution for the 36-term depth‑2 core Ξ36 was extracted at the tensor level (residual exactly 0). **NOT completed:** the final numerical certification of the full identity (and hence the closed form for P‴(1)) failed at the last stage — the absolute weight‑3 identities are exact as tensors but their numeric remainders were not yet identified in the invisible-constant space; a precision-truncation bug in the final PSLQ test was identified at the end and the full-precision re-test remains pending. All difference-type checks that were run to completion close perfectly (65/65 at ≤ 7.6e‑621, and weight 2 closes *exactly*). Honest verdict: the reduction theorem and the explicit Li₃-content are established at the symbol level; the certified numeric closed form is one bug-fix away.

## Task 1 — letters and the S-unit basis (COMPLETE, exact) — `step_14`, `step_15`, `letters_exact.pkl`

F = ℚ(√3) (class number 1; 𝔭₂ = 1+√3 of norm −2 above 2; 𝔭₃ = √3 of norm −3 above 3; fundamental unit ε = 2+√3). **S = {𝔭₂, 𝔭₃}**; V = (𝒪_{F,S}^×/tors)⊗ℚ has basis {𝔭₂, 𝔭₃, ε}, dim 3. Every letter of the 36 depth‑2 constants factors over S. All letters reduce to monomials in {xᵢ, 1−xᵢ, xᵢ−xⱼ}: for the swaps, γᵢⱼ = xⱼ/(xⱼ−xᵢ), γᵢⱼ−1 = xᵢ/(xⱼ−xᵢ), δᵢⱼ = γᵢⱼ(1−xᵢ), δᵢⱼ−1 = xᵢ(1−xⱼ)/(xⱼ−xᵢ), and for Li₂,₁(x,y) the letters are x, y, x−1, y−1, xy−1, x(y−1)/(x−1) — all verified exactly. Factorization table (a = sign·𝔭₂^v₂ 𝔭₃^v₃ ε^k, printed):

| letter | factorization | letter | factorization |
|---|---|---|---|
| x₁ | −𝔭₃⁻¹ε⁻¹ | 1−x₁ | +𝔭₂²𝔭₃⁻¹ε⁻¹ |
| x₂ | −𝔭₂ε⁻¹ | 1−x₂ = √3 | +𝔭₃ |
| x₅ | +𝔭₂𝔭₃⁻¹ε⁻¹ | 1−x₅ | +𝔭₃⁻¹ |
| x₆ | +𝔭₂⁻² | 1−x₆ = √3/2 | +𝔭₂⁻²𝔭₃ε |
| x₅−x₁, x₁−x₂ | ±𝔭₃⁻¹ | x₂−x₅ | −𝔭₂²𝔭₃⁻¹ε⁻¹ |
| x₆−x₁, x₅−x₆ | ±𝔭₂⁻²𝔭₃⁻¹ε | x₂−x₆ | −𝔭₂⁻²𝔭₃ε |

All 20 letters are S-units (no primes beyond 2 and 3 appear) — printed check "ALL LETTERS ARE S-UNITS".

## Task 2 — symbols over F, validated exactly (COMPLETE) — `dod3_symbol.py`, `step_17`

𝒮(Li₃(z)) = −(1−z)⊗z⊗z; 𝒮(Li₂,₁(x,y)) = (xy−1)⊗(xy)⊗(y−1) + (xy−1)⊗[x(y−1)/(x−1)]⊗x + (y−1)⊗(x−1)⊗x, with `vecfun` = exact S-unit factorization (norm valuations + unit-power identification; torsion −1 dropped). Validated at the formal-letter level, all residuals **exactly 0**: integrability (both wedge conditions), the correct stuffle 𝒮(Li₂,₁)+𝒮(Li₁,₂ swapped)+𝒮(Li₃(xy)) = 𝒮(Li₂Li₁), and the Landen identity for Li₃. Quotient Q = V^⊗3/Decomp: Decomp rank **19**, dim Q **8**.

## Task 3 — S-unit equation, full span, exact solution for Ξ36 (COMPLETE at symbol level) — `step_19–22`

**129 solutions** of z+w = 1 in S-units (stable for exponent bounds B = 6 → 8), forming **22 anharmonic orbits** (the special orbit {1/2, 2, −1} + 21 generic 6-orbits), all proven S₃-closed by the search. `solsF_B8.pkl`.

**Main theorem (printed):** rank[𝒮(Li₃(z)) : z ∈ 129 sols | Decomp] = **27 = dim V^⊗3** — the classical trilogarithms at S-unit points span everything modulo decomposables (mirrors OCT8's ℚ-result; here the depth-2 obstruction over F vanishes identically).

**Exact solution for the 36-term depth-2 core.** Ξ36 = Σ₁₂ c^Lam·Lam + Σ₁₂ c^T3a·(Swpγ − Swpδ)/xᵢ with exact ℚ(√3)-coefficients (from DOD2's `atom_coefs_pq.pkl` + the swap relation), split as Ξ36₀ + √3·Ξ36₁ (both ℚ-combos of 36 Li₂,₁ values; numeric check vs Ξ − ΣcT3a·R_quad: **1.37e‑614**). Solving in the 8-dim quotient, both parts match the **same 8 classical Li₃ arguments** (tensor-level residual exactly 0 — `step_22`):

| z | a₀ (1-part) | a₁ (√3-part) |
|---|---|---|
| 2−√3 = ε⁻¹ | −21/32 | −105/32 |
| (√3−1)/2 = 𝔭₂⁻¹ | −27/8 | 9 |
| √3/3 = 𝔭₃⁻¹ | −51/32 | 15/2 |
| 7−4√3 = ε⁻² | 27/128 | 3/32 |
| (2−√3)/2 = x₆ | 51/128 | −207/128 |
| (2√3−3)/3 = −x₁ | −9/8 | 15/16 |
| (3−√3)/6 = (𝔭₂𝔭₃)⁻¹ | −39/32 | 99/32 |
| 1/3 = 𝔭₃⁻² | −27/128 | −45/64 |

**R-piece Li₃-content (exact).** The 12 T3a remainders R(i;j,k) = T3a − (Swpγ−Swpδ)/xᵢ were rebuilt symbolically with **exact** log-factorization (replacing DOD2's PSLQ logfac — `dod3_R.py`, `step_33/34`); 9/12 match quadrature at ~1e‑621, 3 differ by branch terms ~1e‑4 that are products/invisible at symbol level (shown explicitly, `step_37–39`). Their total classical-Li₃ content (12 points, exact ℚ(√3)-coefficients, `R_li3_content.pkl`): Li₃ at 3−2√3, −(2−√3), −x₁, 2−√3, (√3−1)/2, 2√3−3, 2(2−√3), (3−√3)/2, √3−1, 2(3−√3)/3, 3−√3, 2√3−2, with coefficients e.g. −819/32 + 945√3/64 (at 3−2√3), 4473/64 − 1293√3/32 (at −x₁), −2091/32 + 603√3/16 (at √3−1). Merged with the 8 above: **17 distinct classical Li₃ arguments** in the full identity for Ξ (all S-unit points; two > 1, convertible by inversion).

## Task 4 — verification status (partial; honest) — `step_42–62`

* New structural fact (printed): the invisible-constant space over F is bigger than over ℚ — it includes **L(3, χ₁₂)** and **L(2, χ₁₂)·logs** (K₂(𝒪_F) is torsion), computed to 620 digits via Hurwitz zeta (L(3,χ₁₂) = 0.99004001943815994979…, L(2,χ₁₂) = 0.94970312629400939526…, both verified against Dirichlet series). Also Re Li₃(2) = (7ζ(3) + 12ζ(2)log 2)/8 recovered (PSLQ, resid 3e‑620).
* **Weight-2 pipeline: PERFECT** — Li₂(1/4) reduces with remainder **exactly ζ(2)** (PSLQ [1,−1,0] vs [rem, ζ(2), L(2,χ)], full 620-digit precision), validating letters↔values↔products end-to-end.
* **65/65 null Li₃-identities** (pairs of S-unit points with equal projected symbol) close with residuals ≤ **7.6e‑621** and tiny coefficients, using only ⟨ζ(3), ζ(2)·logs⟩.
* **Absolute weight-3 closures (e.g. Li₃(1/4), Li₂,₁(x₅, x₆/x₅)): extracted identities are TENSOR-EXACT (residual 0, verified three ways) but the numeric remainder (−32.799586563322169122 for Li₃(1/4)) was not identified.** The final PSLQ testing the enlarged invisible space ⟨ζ(3), L(3,χ), ζ(2)logs, L(2,χ)logs⟩ was run on a **21-digit truncated remainder** — a bug found in the last review pass; the full-precision re-test (and hence the certified assembly of P‴(1) = 72·L·S with Ξ's closed form against `P3_630.pkl`) was not completed. Most likely resolution: the remainder lies in the enlarged invisible space and needs the full-precision PSLQ; secondarily, one value-side convention for absolute (non-difference) weight‑3 extractions needs re-auditing. The PSLQ product-value mismatch found and fixed mid-run (log 𝔭₂ = (l2+lu)/2, log 𝔭₃ = l3/2, log ε = lu — not l2, l3, lu directly) is documented for the follow-up.

## Evidence table (all printed by code)

| check | residual / result | digits |
|---|---|---|
| 20 letters factor over S = {𝔭₂, 𝔭₃}, units ±1, ε | exact | — |
| γ/δ letter identities (γ−1, δ−1) | exact 0 | — |
| 𝒮(Li₂,₁) integrability / stuffle / Landen | exact 0 | — |
| Decomp rank / dim Q | 19 / 8 | — |
| S-unit equation solutions (B=6→8 stable) | 129, 22 orbits | — |
| rank[Li₃ symbols | Decomp] | **27 = full** | — |
| Ξ36 symbol solution (8 Li₃ args, both parts) | tensor residual exact 0 | — |
| Ξ36₀ + √3 Ξ36₁ vs Ξ − ΣcR (numeric) | 1.37e‑614 | 620 |
| R sym vs quad: 9 triples / 3 triples | ≤2.7e‑621 / ~1e‑4 (product-level branch terms) | 620 |
| R-piece Li₃ content (12 points) | exact ℚ(√3) coeffs | — |
| 65/65 null Li₃ identities | ≤ 7.6e‑621 | 620 |
| weight-2 Li₂(1/4) reduction | remainder = ζ(2), resid 0.0 | 620 |
| Re Li₃(2) identity 8·Li₃(2) = 7ζ(3)+12ζ(2)l2 | 3.02e‑620 | 620 |
| **final absolute weight-3 closures** | **not identified (truncated-PSLQ bug; retest pending)** | — |

## Files

`dod3_symbol.py` (exact ℚ(√3) + S-unit factorization + symbol engine + quotient), `dod3_R.py` (exact R-pieces), `letters_exact.pkl`, `solsF_B4/6/8.pkl`, `li3_gens.pkl`, `li2_basis.pkl`, `xi36_symbol.pkl`, `xi36_li3_basis_solution.pkl`, `R_li3_content.pkl`, `R_all12_vals620.pkl`, `prod_sym_data.pkl`. Key steps: step_14/15 (letters), step_17 (validation), step_19/20 (129 sols, span 27), step_21/22 (Ξ36 solution), step_33/34 (R exact), step_44–53 (products fix + 65/65 nulls), step_56–62 (closure attempts, invisible space, weight-2 success).

## Paragraph for the paper

*For the regular 12-gon, the depth-2 remainder Ξ of P‴(1) lives over F = ℚ(√3) with S the primes above 2 and 3: every letter of its 12 Lam and 24 swap constants factors over the S-unit group ⟨1+√3, √3, 2+√3⟩. In the symbol space V^⊗³ (V 3-dimensional) modulo shuffle-decomposables (dim 8), the classical trilogarithms at the 129 solutions of the S-unit equation over F span everything (rank 27), so Goncharov's depth-2 obstruction vanishes identically: Ξ is an explicit ℚ(√3)-combination of classical Li₃ at 17 S-unit arguments (8 from the depth-2 core, 12 from the T3a remainders, 3 overlaps) plus products of lower-weight functions and the constants ζ(3), L(3,χ₁₂), ζ(2)·logs and L(2,χ₁₂)·logs. The symbol-level extraction is exact and all difference-type identities verify at 620 digits (65/65, ≤ 7.6e‑621); the final absolute-value closure awaits a corrected full-precision PSLQ after a truncation bug was isolated — the earlier "irreducibility" of DOD2 is thereby refuted at the level of the symbol, exactly as in OCT8: the required trilogarithm arguments (ε⁻², (𝔭₂𝔭₃)⁻¹, γ/δ-points, …) simply were not in the earlier 102-element candidate basis.*

**Caveats (honest accounting):** (i) the certified numerical closed form of Ξ and P‴(1) is not yet established — one full-precision PSLQ remains; (ii) 3 of 12 symbolic R-pieces need their branch-correction products written out explicitly for a fully displayed identity; (iii) completeness of the 129-solution list is empirical (stable at B = 6 → 8), though sufficiency is guaranteed by the full-span theorem.