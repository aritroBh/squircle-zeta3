# DEC2 complete — the regular decagon's P'''(1) closed form is **PROVED (computer-assisted)** end to end

All four planner items are done. **Nothing stays numerical.** `verify_dec2.py` re-runs the entire chain in **169 s (< 5 min): ALL PARTS PASSED (P1–P6)**.

## The result
**P'''(1) = 60·L·S**, L = √((25−11√5)/10) (edge length; L = ρ(3√5−5)/20 in ℚ(ρ), ρ = √(10−2√5)), with S the 15-term ℚ(√5)-combination (rational term exactly 0):

| constant | p | q | | constant | p | q |
|---|---|---|---|---|---|---|
| Li₃(1/φ) | −75/32 | 81/64 | | log 2 | 1 | 0 |
| Li₃(−1/φ) | −35/32 | 53/64 | | log φ·ζ(2) | 3/16 | −13/320 |
| ζ(3) | 33/64 | −19/80 | | log²φ | −73/8 | 29/8 |
| ζ(2)=π²/6 | −1/16 | −9/80 | | log 2·ζ(2) | 1/32 | 9/160 |
| log φ | 5/2 | −1/2 | | log 2·log φ | −5 | 1 |
| log²2 | −1 | 0 | | log³φ | 1099/96 | −1991/384 |
| log 2·log²φ | 81/16 | −33/16 | | log²2·log φ | 5/4 | −1/4 |
| log³2 | 1/6 | 0 | | | | |

Exact match (0 mismatches, term by term) with DEC1's `final_table_pq.pkl` and with the paper's prose (15 terms over ℚ(√5), no depth 2, rational term 0).

## What was proved and how
**Task 1 — exact integrand (PROVED).** Geometry derived from scratch over ℚ(ρ): normals (cos(iπ/5), sin(iπ/5)); edge V₀ = (−h, y₀), V₁ = (−h, −y₀) with h = (√5−1)/4 = 1/(2φ), y₀ = (√5−2)/ρ; kink conditions exact; Σaᵢ(t) = 1 exact; term values a₂ = t/(2φ²), a₃ = (1−t)/(2φ²), a₁ = (1+t/φ)/(2φ²), a₄ = (1−t/φ²)/(2φ), a₅ = 1/(2φ), all coefficients in ℚ(√5). Length element from the proved radial map x = F^{−1/(1+d)}y: |x'|² = F^{−2/(1+d)}(|y'|² − 2(y·y')Φ + |y|²Φ²). Exact series engine over ℚ(√5) gives E₃n = [d³]√(|x'|²/L²): 55 log-monomials; **0 mismatches vs DEC1's reconstructed E₃n (440 rational coefficients)**; independent direct d-expansion at t = 0.3: |diff| = 1.03e-80. Only radical outside ℚ(√5): the overall L.

**Task 2 — exact integrals (PROVED).**
- Local (40): antiderivative certificates F′ = f (formal identities) + exact endpoints; 0 mismatches vs OCT5's Jloc.
- Single-ratio J-table (40×2 sides): exact symbolic summation (esum machinery itself justified: 3 generating functions, 3 atom identities by x·d/dx induction, w^c expansions — residuals exactly 0); 0 mismatches vs K1d; quadrature worst 2.1e-99.
- D1–D10: D1 = Λ/x₁, D2 = swap/x₄ by absolutely-convergent series rearrangements (term equality exact); D3, D4 by exact t→1−t reflections (5 algebraic identities = 0); D5–D8 by explicit antiderivatives (canonical-letter certificates, all coefficients exactly 0) + exact endpoints; D9, D10 by IBP from D7, D8 (residual 0, boundaries ±l_φ³/2 exact). All match dec_engine symbolically; quadrature worst 2.5e-111.
- 20 mixed integrals: 20 exact IBP integrand identities (residuals 0; M₂ uses the regularized primitive (t^{n+1}−1)) + exact partial fractions; 0 mismatches vs engine; quadrature worst 4.4e-110.

**Task 3 — exact assembly (PROVED).** 261 monomials (0 mismatches vs DEC1's mons.pkl) → S in the named basis → proved reductions: golden Li₂'s (Euler+Landen+duplication, certified by differentiation+pins); Li₃(1/φ²) = 4Li₃(1/φ)+4Li₃(−1/φ) (duplication, exact series); Lewin via the certified 3-term Landen identity Li₃(x)+Li₃(1−x)+Li₃(x/(x−1)) = ζ(3)+ζ(2)log(1−x)−½log²(1−x)log x+⅙log³(1−x) (derivative residual 0); A-atoms (DEC1/OCT6 orbit system — the L1/L2/DIFF derivative residuals re-verified here: exactly 0); Λ and swap (PRV1 Theorem D — verify_full.py re-run: ALL CHECKS PASSED, 11 s). Result: the 15-term table, **0 mismatches vs DEC1**; **local ζ(3) = C₅ exactly** (60L(5/32−√5/8) − C₅ = 0 in ℚ(ρ); C₅ = ρ(105√5−255)/32 matches Theorem C_m at m = 5); **60·L·S vs the 676-digit gold: |diff| = 3.9e-662** (700-dps evaluation).

## Evidence table (all numbers printed by code)
| check | label | residual |
|---|---|---|
| E₃n symbolic vs DEC1 (55 keys, 440 coeffs) | PROVED | 0 mismatches |
| E₃n vs direct d-expansion (80 dps) | check | 1.03e-80 |
| 40 local integrals (certificates + endpoints) vs OCT5 Jloc | PROVED | 0 mismatches |
| J-table (esum) vs K1d, 40 keys | PROVED | 0 mismatches |
| J-table at x₁, x₄ vs quadrature (100 dps, 80 values) | check | 2.1e-99 |
| D1–D10 vs dec_engine | PROVED | 0 (symbolic) |
| D1–D10 vs quadrature (110 dps) | check | 2.5e-111 |
| 20 mixed integrals: IBP identities + engine match | PROVED | 0 (residuals/symbolic) |
| 20 mixed vs quadrature (110 dps) | check | 4.4e-110 |
| Classical identities (Euler/Landen/inversion/Li₃-Landen/A-atom L1,L2,DIFF) | PROVED | 0 (derivative residuals) |
| Golden Li₂ ×3, Li₃(1/φ²) dup + Lewin, sum relation (120 dps) | PROVED | ≤1.5e-121 |
| A-atom reductions ×2 (DEC1/OCT6 proofs, 120 dps) | PROVED (cited) | ≤4.9e-122 |
| Λ, swap reductions (PRV1 Theorem D; verify_full.py re-run) | PROVED (cited) | ALL CHECKS PASSED |
| 15-term table vs DEC1 final_table_pq.pkl | PROVED | 0 mismatches |
| local ζ(3): 60L·c − C₅ | PROVED | 0 exactly |
| **60·L·S vs 676-digit gold (700-dps eval)** | check | **3.9e-662** |
| verify_dec2.py (self-contained, end to end) | — | ALL PARTS PASSED, 169 s |

## Deliverables (work dir `new_results/DEC2_decagon_exact_assembly/`)
`verify_dec2.py` (25.2 KB, runs 169 s), `result.md`, `progress.log`, `esum_engine.py`, and pickles `geom_exact.pkl`, `E3n_sym_pq.pkl`, `Jloc_dec2.pkl`, `Jloc_antiderivatives.pkl`, `Jtab_dec2.pkl`, `D_all_proved.pkl`, `MIX20_proved.pkl`, `mons_dec2.pkl`, `S_assembly_dec2.pkl`, `final_table_dec2.pkl`; session trail in `step_01.py…step_99.py` (planner transcript).

**Status: the decagon row of Table II can be upgraded from "Numerical (reductions proved)" to PROVED (computer-assisted): exact integrand, exact integrals, exact assembly — verified end to end by `verify_dec2.py` (ALL PARTS PASSED), with the 676-digit quadrature agreeing to 3.9e-662.**