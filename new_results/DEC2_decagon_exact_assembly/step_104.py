import os
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
md = """# DEC2 — the regular decagon's P'''(1): exact integrand and exact assembly — **PROVED (computer-assisted)**

**Bottom line.** Every step of the decagon's closed form is now exact:
  P'''(1) = 60 L S,  L = sqrt((25-11 sqrt5)/10)  (the edge length; L = rho(3 sqrt5 - 5)/20 in Q(rho),
  rho = sqrt(10-2 sqrt5)),
with S the 15-term Q(sqrt5)-combination of {l2=log 2, lf=log phi, z2=pi^2/6, z3=zeta(3),
Li3(1/phi), Li3(-1/phi)} tabulated below (rational term exactly 0). The integrand E3 is derived
symbolically (0 mismatches vs DEC1's reconstruction, 440 rational coefficients); every integral
is proved (antiderivative certificates, exact series, reflections, IBP); every reduction is
proved (golden Li2/Li3 lattice by differentiation+pins; A-atoms by DEC1/OCT6's proved function
identities; the two depth-2 constants by PRV1's proved theorems); the assembly is exact and
agrees with the 676-digit quadrature to 3.9e-662. Nothing remains numerical.
`verify_dec2.py` re-checks everything end to end in 169 s: **ALL PARTS PASSED**.

## The 15-term closed form (matches DEC1 final_table_pq.pkl and the paper exactly; (p,q) = p + q sqrt5)

| constant | p | q |
|---|---|---|
| Li3(1/phi)        | -75/32 | 81/64 |
| Li3(-1/phi)       | -35/32 | 53/64 |
| z3                | 33/64 | -19/80 |
| z2                | -1/16 | -9/80 |
| lf                | 5/2 | -1/2 |
| l2                | 1 | 0 |
| lf*z2             | 3/16 | -13/320 |
| lf^2              | -73/8 | 29/8 |
| l2*z2             | 1/32 | 9/160 |
| l2*lf             | -5 | 1 |
| l2^2              | -1 | 0 |
| lf^3              | 1099/96 | -1991/384 |
| l2*lf^2           | 81/16 | -33/16 |
| l2^2*lf           | 5/4 | -1/4 |
| l2^3              | 1/6 | 0 |

(rational part: exactly 0.) P'''(1) = -4.776200180882238416632...

## Task 1 — exact integrand (PROVED)
- Geometry from scratch over Q(rho) (rho = sqrt(10-2 sqrt5), degree 4, contains sqrt5): normals
  (cos(i pi/5), sin(i pi/5)); edge V0 = (-h, y0), V1 = (-h, -y0), h = (sqrt5-1)/4 = 1/(2 phi),
  y0 = (sqrt5-2)/rho; kink conditions exact; sum_i a_i(t) = 1 exact. Term values:
  a2 = t/(2 phi^2), a3 = (1-t)/(2 phi^2), a1 = (1+t/phi)/(2 phi^2), a4 = (1-t/phi^2)/(2 phi),
  a5 = 1/(2 phi) — all coefficients in Q(sqrt5) (normal-form check). L^2 = (25-11 sqrt5)/10 exact.
- Length element (proved geometry): radial map x(t) = lambda y(t), lambda = F^{-1/(1+d)},
  F = sum a_i^{1+d}; |x'|^2 = F^{-2/(1+d)} (|y'|^2 - 2 (y.y') Phi + |y|^2 Phi^2), Phi = S/F.
- Exact series engine over Q(sqrt5): E3n = [d^3] sqrt(|x'|^2/L^2): 55 log-monomials, t-degree <= 3.
  Monomial-by-monomial vs DEC1's E3n_exact_pq.pkl: 55/55 keys, **0 mismatches**.
- Independent direct numerical d-expansion at t = 0.3 (80 dps): |diff| = 1.03e-80.
- Only radical outside Q(sqrt5): the overall edge length L (the analog of the octagon's u).

## Task 2 — exact integrals (PROVED)
- **Local (40 keys)**: antiderivative certificates F' = f verified as formal identities + exact
  endpoint evaluation; values match OCT5's Jloc symbolically (0 mismatches).
- **Single-ratio J-table (40 keys)**: exact symbolic summation (esum machinery justified: 3
  generating functions, 3 atom identities by x d/dx induction, w^c expansions c=2,3 by
  differentiation — all residuals exactly 0; moment identities = Beta/psi formulas, also certified
  at l=1..4 by the local table); rebuilt table matches K1d's Jtab_t.pkl symbolically (0 mismatches);
  validated at x1 = -1/phi, x4 = 1/phi^2 vs quadrature: worst 2.1e-99.
- **D1..D10 master constants**:
  D1 = Lam/x1, D2 = swap/x4: absolutely convergent double-series rearrangements (term-by-term
  equality exact; |x1| < 1, x4 < 1);
  D3, D4: t -> 1-t reflections of D2, D1 (5 exact algebraic substitution identities);
  D5..D8: explicit antiderivatives (F5..F8, verified as formal identities by the canonical-letter
  method: every letter coefficient exactly 0) + exact endpoint evaluation;
  D9, D10: IBP from D7, D8 (integrand identities residual 0; boundaries -lf^3/2, +lf^3/2 exact).
  All match DEC1's dec_engine symbolically (diff 0); all vs quadrature at 110 digits: worst 2.5e-111.
- **20 mixed integrals** (5 types x n = 0..3): each reduced by an exact IBP identity
  (20 integrand identities, residuals 0; M2 uses the regularized primitive (t^{n+1}-1) with both
  boundary limits exactly 0) + exact partial-fraction identities; values match dec_engine
  symbolically (0 mismatches); vs quadrature at 110 digits: worst 4.4e-110.

## Task 3 — exact assembly (PROVED)
- mons (261 log-monomials in the reduced basis) rebuilt from my E3n: 0 mismatches vs DEC1's mons.pkl.
- Assembly S = sum coef * INT over the 261 monomials (named basis {l2, lf, z2, z3, Li2/Li3/A at
  x1, x4, 1/phi, Lam, swap}), then the PROVED reductions:
  golden Li2: Li2(1/phi^2) = 2 z2/5 - lf^2, Li2(1/phi) = 3 z2/5 - lf^2, Li2(-1/phi) = -2 z2/5 + lf^2/2
    (linear algebra from Euler + Landen + duplication for Li2 — all three certified by
    differentiation + exact pins);
  Li3(1/phi^2) = 4 Li3(1/phi) + 4 Li3(-1/phi)  (duplication, exact series);
  [consequence, Lewin: Li3(1/phi^2) = 4(z3 - lf z2)/5 + 2 lf^3/3, from the certified 3-term
    Landen identity Li3(x)+Li3(1-x)+Li3(x/(x-1)) = z3 + z2 log(1-x) - (1/2) log^2(1-x) log x
    + (1/6) log^3(1-x) — derivative residual exactly 0, pin at x -> 0+ — plus duplication];
  A-atoms A(-1/phi) = 4 z3/5 - 2 lf z2/5 + 2 Li3(-1/phi),
          A(1/phi^2) = z3 - 3 lf z2/5 + 4 Li3(-1/phi) + 3 Li3(1/phi)
    (DEC1's orbit system of the three function identities L1/L2/DIFF — proved in OCT6 by
    differentiation + pins, and re-verified here: all three derivative residuals exactly 0);
  depth-2: Lam = Li21(-1/phi,-1/phi) = 2 z3/5 + 4 lf z2/5 - lf^3 + 3 Li3(-1/phi) + Li3(1/phi),
           swap = Li21(1/phi^2,-phi) = 3 z3/5 - 7 lf z2/5 + 2 lf^3 + Li3(-1/phi) + Li3(1/phi)
    (PRV1 Theorem D — proved there; verify_full.py re-run here: ALL CHECKS PASSED, 11 s).
  All reduction statements additionally validated at 120 digits: worst 1.5e-121.
- Result: the 15-term table above — **0 mismatches** vs DEC1's final_table_pq.pkl (term-by-term,
  exact); rational term exactly 0; agrees with the paper's prose (15 terms over Q(sqrt5), no depth 2).
- **Local zeta(3) = C5 exactly**: the local-only part carries z3-coefficient 5/32 - sqrt5/8 per
  edge-factor, and 60 L (5/32 - sqrt5/8) - C5 = 0 exactly in Q(rho), with
  C5 = rho (105 sqrt5 - 255)/32 = the C_m of Theorem "Regular 2m-gons" at m = 5 (exact match).
- **Final numeric check**: 60 L S vs DEC1's 676-digit tanh-sinh gold P'''(1):
  |diff| = 3.9e-662 (evaluated at 700 dps); in verify_dec2.py at 150 dps: 8.0e-148.

## Status by step (labels)
- Task 1 (exact integrand): **PROVED**.
- Task 2 (exact integrals): **PROVED**.
- Task 3 (exact assembly + 676-digit agreement): **PROVED** (the 676-digit value itself remains
  DEC1's numerical gold standard, used only as evidence; the closed form is proved independently).
- Nothing stays numerical.

## Files
step_01..step_15 (this session's trail in the planner transcript), verify_dec2.py (self-contained,
169 s, ALL PARTS PASSED), geom_exact.pkl, E3n_sym_pq.pkl, Jloc_dec2.pkl, Jloc_antiderivatives.pkl,
Jtab_dec2.pkl, D_all_proved.pkl, MIX20_proved.pkl, mons_dec2.pkl, S_assembly_dec2.pkl,
final_table_dec2.pkl, esum_engine.py, progress.log, result.md.
"""
with open("result.md","w") as f:
    f.write(md)
with open("progress.log","a") as f:
    f.write("step_16: verify_dec2.py written and run (169s): ALL PARTS PASSED (P1..P6). result.md finalized: decagon P'''(1) closed form PROVED end to end. DONE.\n")
print("result.md finalized,", len(md), "chars")
