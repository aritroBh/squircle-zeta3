# OCT6 — The octagon's P'''(1) in closed form

## Headline result

**P'''(1) = 48·u·S**, u = √(2+√2), where S is an explicit ℚ(√2)-combination of the finite named set
{ l2 = log 2, lu = log(1+√2), z2 = π²/6, z3 = ζ(3), and Li₂, Li₃ at the six nontrivial real reps
1−√2, 1−√2/2, √2−1, 2−√2, √2/2, −√2/2 } — **no A-atoms and no depth-2 constants**: Λ = Li₂₁(1−√2, −√2/2)
and its swap both reduce to single Li₃'s. Verified **3.6e-412** against an independent 420-dps quadrature
of the entire integrand (the stored OCT5 P3_720 value is only ~203 digits accurate — its quad was
under-resolved; we match all of its accurate digits, 2.8e-200).

**The 35-term closed form** (P'''(1)/(48u); rational part is exactly 0):

```
w1 l2: 5-7√2/2 ; lu: 3-2√2
w2 L2d: -3/4+√2/2 ; L2k: 3/4-√2/2 ; L2x1: -7/4+5√2/4 ; L2x2: 7/4-5√2/4 ;
   l2·lu: -77/8+27√2/4 ; l2²: -389/32+275√2/32 ; lu²: -3+2√2 ; z2: -3/4+√2/2
w3 L3d: 3/8-√2/4 ; L3k: -3/8+√2/4 ; L3mr2h: -9/2+3√2 ; L3r2h: -155/32+103√2/32 ;
   L3x1: 33/32-21√2/32 ; L3x2: -79/32+59√2/32 ; l2·L2d: -17/16+49√2/64 ;
   l2·L2k: 11/32-15√2/64 ; l2·L2mr2h: 41/32-7√2/8 ; l2·L2r2h: -41/32+7√2/8 ;
   l2·L2x1: -5/32+3√2/32 ; l2·L2x2: -41/16+57√2/32 ; l2·lu²: 63/32-87√2/64 ;
   l2·z2: 101/32-71√2/32 ; l2²·lu: 443/64-631√2/128 ; l2³: 8635/768-3053√2/384 ;
   lu·L2d: 3/8-√2/4 ; lu·L2k: -3/8+√2/4 ; lu·L2r2h: -5/32+√2/8 ;
   lu·L2x1: 63/32-45√2/32 ; lu·L2x2: 1/16-√2/32 ; lu·z2: 17/32-3√2/8 ;
   lu³: 1/2-√2/3 ; z3: 29/32-39√2/64
```
(x1=1−√2, x2=1−√2/2, k=√2−1, d=2−√2, r2h=√2/2, mr2h=−√2/2; table_Stot.pkl)

## How it was proved

1. **83 monomials** (step_01): every E3 coefficient is u·(p+q√2)·tⁿ exactly (1/√(10−7√2) = u(2+√2)/2; reconstruction error 9.5e-33). Classes: 34 local, 21+21 single-ratio, **7 mixed types**.
2. **Mixed part exactly** (Task 2, steps_02–04): IBP/partial fractions reduce all mixed integrals to the K1d J-table (verified at x1 and x2), Jloc, and ten master constants D1–D10, all of which close: D1=−Λ/κ1, D2=−swap/κ2 (exact series identity ∫₀¹ X log(1−bt)/(1−at)dt = Li₂₁(a,b/a)/a); D3,D4 → Λ, swap, Li₂(x1), Li₂(x2) (via t↦1−t); D5–D8 → Li₂/Li₃ at {√2−1, 2−√2} + log products (exact antiderivatives); D9,D10 → D7,D8 + l2³ (IBP). Each D verified by quadrature (≤2.6e-131); all 20 reduced mixed integrals verified 1.98e-103. Mixed = u·S_mix = 0.00169887193388033285629748085803… (fresh 760-dps quad: 1.1e-576).
3. **A-atoms** (Task 1, steps_10,19): A(x) = Nielsen S₁,₂. Three function identities, each derived by differentiation (derivatives collapse by dilog Landen/Euler/inversion) and point-pinned — Landen-sum L2 (745 digits ×3 pts), inversion L1 (400 digits ×3 pts), and the new difference identity
   **A(x)−A(1−x) = 2Li₃(x)−2Li₃(1−x)+log(1−x)Li₂(1−x)−log x·Li₂(x)−½log²x·log(1−x)+½log²(1−x)·log x** (400 digits ×5 pts) —
   form a solvable 4×4 orbit system: A(1−√2) = [R_L2(x2)−DIFF(x2)−R_L2(√2/2)+R_L1(κ1)]/2, A(1−√2/2) = R_L2(x2)−A(x1). These **prove** the LLL reductions R1, R2 to the 8-rep Li₃ basis (match direct series 7e-746).
4. **Λ and swap** (R3, R4): LLL-discovered, verified 745 digits; R3 independently cross-validated by the stuffle Li₂(x1)Li₁(−√2/2) = Λ + Li₁₂(−√2/2, x1) + Li₃(x2) with Li₁₂ from a fresh double series (6.7e-744).
5. **Assembly & independent basis** (Task 3, 4): full assembly → 35-term table (412-digit verified). Relation lattice (LLL among named constants only): 5 weight-2 relations (all Li₂'s → one Li₂ + logs/π²) and 13 weight-3 relations (two Landen triples T1{κ1,d,mr2h}, T2{x2,r2h,x1}, the duplication L3r2h+L3mr2h = 7z3/32−z2·l2/8+l2³/24, plus 10 induced Li₂·log relations) — every one a classical identity independently re-verified at ≥745 digits. Independent basis: 8 constants (w≤2) + 12 (w3).

## Requested reports

- **Coefficient of ζ(3).** Basis-free: the **local** part (Jloc only, where ζ(3) enters only explicitly) has ζ(3)-coefficient 48·(9/8−13√2/16) = **54−39√2 — exactly C4/u**. In the named presentation, total ζ(3)-coefficient of P'''(1)/u is **87/2 − 117√2/4**, so the **non-local ζ(3) is N_e = (39√2−42)/4 per edge** (×u). The local/non-local split is canonical (the local part carries no Li₃'s); the absolute coefficient beyond that is **basis-dependent** — e.g. after eliminating via the z3-bearing Li₃ relations it reads −55891/512 + 39501√2/512 per 48u.
- **Depth-2 survival: none.** Λ and swap (hence all Li₂₁/Li₁₂/Li₁₁ at the octagon points) reduce to single Li₃'s. OCT5's LLL obstruction is explained: the needed rep-Li₃ directions (L3k, L3d, …) lie outside the old 20-element basis (its Li3a/b/c are all in the rep span — checked 745 digits), so Λ/swap only *appeared* to be new constants.
- **Difference from C4 = (54−39√2)u.** P'''(1)/u − C4/u = 48·S − (54−39√2) (table_diffC4.pkl, 1.7e-93): a full weight-≥1 combination (not a bare ζ(3)-term); its ζ(3)-component is N_e = (39√2−42)/4. Equivalently P'''(1) = C4·ζ(3) + 48u·(S_s1+S_s2+S_mix) exactly.

## Evidence table

| check | residual |
|---|---|
| 83 monomial coefficients = u·(p+q√2) | 9.5e-33 |
| Jtab (40 keys) at x1, x2; Jloc (40) | 4.8e-52, 4.8e-51 |
| D1–D10 closed forms vs quadrature (130 dps) | ≤ 2.6e-131 |
| 20 reduced mixed integrals vs quadrature | 1.98e-103 |
| R1, R2 (A-atoms) — **proved**, vs direct series | 7.3e-746, 1.1e-746 |
| R3 (Λ), R4 (swap) vs 745-digit eval | 6.6e-744, 5.3e-744 |
| R3 via stuffle + independent Li₁₂ series | 6.7e-744 |
| 5 Li₂ + 13 Li₃ lattice relations (classical) | ≤ 2.2e-746 |
| mixed part vs fresh 760-dps quad | 1.14e-576 |
| **P'''(1) exact vs 420-dps quad of whole E3** | **2.1e-411** |
| P'''(1) vs OCT5 stored P3_720 (~203 good digits) | 2.8e-200 (all of them) |

## Files (work dir OCT6_octagon_closed)

`result.md` (full write-up), `progress.log`; producing steps: **step_01** (monomials → mons_simple.pkl), **step_02** (D1–D10, Dconsts.pkl), **mixed_engine.py** + **step_03** (exact IBP engine, 20 mixed integrals), **step_04** (assembly, S_assembly.pkl), **step_05/07** (quad checks, P3_exact.pkl, 745-digit P'''(1)), **step_06** (Jtab/Jloc table verification), **step_10/11** (reductions.pkl, 745-digit re-verification), **step_14** (relation lattice, indep.pkl), **step_15** (S_indep.pkl, final_tables.pkl), **step_16/17** (classical relations, Li3a/b/c span, table_diffC4.pkl), **step_19** (analytic proof of R1/R2), **step_20** (final end-to-end 3.6e-412). OCT5 inputs preserved in `oct5/`; K1d J-table read-only.