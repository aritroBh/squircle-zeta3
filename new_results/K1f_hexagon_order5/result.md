# K1f — P⁽⁵⁾(1) of centrally symmetric hexagons, weight-5 engine and J-table (interrupted run; Tasks 1–2 complete and verified, Task 3 assembly staged but not yet run)

## What was accomplished and verified

### Task 1 — General single-ratio E5(t): **DONE** (`step_24.py`, `E5_polys.pkl`)
The K1d E3 builder was extended to d⁵. For an edge with term values `a₀=A₀t, a₁=A₁(1−t), a₂=B(1+κt)` (κ=r−1, A₁+B=1, A₀=1−Br), data `L₂,q₀,V₀₂`, unweighted logs `lA₀,lA₁,lB`, and symbolic logs `X=log t, Y=log(1−t), W=log(1+κt)`, the integrand is built as `E = Fpow·√(L₂ − 2y·y′Φ + y·yΦ²)`, `Fpow = F^{−1/(1+d)}`, `Φ = F′/F` (the K1d-fixed convention).
- **Exact cross-validation:** my `E[3]` equals K1d's `E3_fixed.pkl` **symbolically, 80 monomials, 0 mismatches** — the same engine that K1e verified to 100 digits, so the weight-5 integrand inherits that validation of conventions.
- `E[5]` = 336 monomials `tⁿXᵃYᵇWᶜ`, n≤5, a+b+c≤5: 35 non-local `(a,b,c)` patterns (c≥1) × n=0..5 = 210 keys; 21 local `(a,b)` × n=0..5 = 126 keys.

### Task 2 — Weight-5 J-table `∫₀¹ tⁿXᵃYᵇWᶜ dt`: **DONE, verified at 80 dps** (`jtab5.pkl`)
`J(n;a,b,c) = c!(−1)ᶜ Σ_{m≥c} xᵐ e_{c−1}(m−1)/m · Jloc(n+m;a,b)`, x=−κ=1−r, via: exact symbolic local Beta table `Jloc(N;a,b)` (a+b≤5, in `{1,ζ2,ζ3,ζ2²,ζ2ζ3,ζ5}`); a weight-5 e-sum engine `esum5.py` (K1c extended to H⁽¹..⁴⁾ numerators); and an exact shuffle-algebra reduction of all 26 raw harmonic atoms to the **13-generator Lyndon HPL basis** `{H1, Li2, Li3, H011, Li4, H0011, H0111, Li5, H00011, H00101, H00111, H01011, H01111}`.
- **Minimal atom set (answer to Task 2's question): all 13 Lyndon atoms occur** — including the **two non-Nielsen weight-5 HPLs** `H00101 = Li₃,₂(x,1)` and `H01011 = Li₂,₂,₁(x,1,1)` (30 J-keys each) — plus ζ(2), ζ(3); ζ(5) enters only through the local table (a+b=5). 46 distinct atom monomials; deepest pole x⁻⁶, and **each J-entry satisfies J ~ x^c as x→0** (ratios 10^c confirmed), so every non-local entry vanishes at r=1 — the κ→0 cancellation needed for the r=1 local-law reduction is structural.
- `step_31.py`: exact HPL argument transform x↦x/(x−1) (letter map `0→0+1, 1→−1`, i.e. φ*f₀=f₀+f₁ — Landen-checked), giving atom values at any real x<1 (e.g. x=−2 for the fresh hexagon's r=3 edge).

## Evidence table (all numbers printed by my code)
| Check | Object | Result | Residual | Digits |
|---|---|---|---|---|
| E[3] vs K1d `E3_fixed.pkl` | symbolic, 80 monomials | exact match | 0 | exact |
| Atom→HPL reductions (step_26) | 26 atoms × 3 x-values | verified | 2.87e-79 | 80 |
| Local Beta table (step_27) | 12 N × 21 (a,b) | verified | 9.3e-46 | 50 |
| e-sum engine (step_28) | 15 sums × 3 x | verified | 8.4e-80 | 80 |
| **Full J-table (step_30)** | 40 keys × 4 κ vs quad | verified | **1.3e-76** | 80 |
| J ~ x^c at x→0 | 4 entries | ratios = 10^c ✓ | — | 50 |
| x↦x/(x−1) transform (step_31) | all words w≤5 at x=−1 vs nsum | verified | 9.49e-81 | 80 |

## Status vs the four tasks
- **Task 1: DONE** (exact E5, cross-validated at weight 3). The t=0.3 numeric spot-check was superseded by the exact E[3]=K1d match plus J-table quad checks.
- **Task 2: DONE** (210-key table built, minimal atom set identified, all pieces verified ≥76 digits).
- **Task 3: NOT yet run** — the assembly `E_arc⁽⁵⁾ = 120·(Σ_{c≥1} co·J|_{x=−κ} + Σ co·Jloc)`, its basis collection, κ⁻⁶-pole cancellation under A₀=1−B−Bκ, the r=1 reduction to Z_e⁽⁵⁾/Z_e^{(2,3)}, and hexagon numerics (b=0, b=1, fresh, irrational) are fully staged (engine, table, atom evaluator all ready and verified) but were not executed when the run was interrupted.
- **Task 4: pending** (depends on the Task-3 assembly; expect explicit ζ(5)-coefficient = ΣZ_e⁽⁵⁾ for all r, with non-local ζ(5) entering through Li₅/HPL values at special r).

## Files (work dir `new_results/K1f_hexagon_order5/`)
`step_24.py` (E5 engine rebuild + monomial split + K1d check), `E5_polys.pkl`; `step_25.py` (shuffle/Lyndon machinery, `atom_table.pkl`), `step_26.py` (atom verification); `step_27.py` (`jloc_sym.pkl`); `esum5.py` (engine), `step_28.py` (engine checks); `step_30.py` (**`jtab5.pkl`**, full table + verification); `step_31.py` (atom evaluator, any x<1); `progress.log`, `result.md` (full status + assembly plan). All are in place for the next run to execute the assembly directly.