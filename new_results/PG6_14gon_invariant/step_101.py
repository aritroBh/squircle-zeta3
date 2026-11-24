summary = r"""# PG6: 14-gon (2,1) invariant — running summary

## Task A — P'''(1) discrepancy: RESOLVED
- CORRECT: **P'''(1) = -5.5080836980534783** (exact-integrand chain 84*L*int E3). Confirmed to ~24 digits
  by 3 methods (high-prec int E3n; Fourier filter on sliver-free chamber; Richardson FD). The polar gold
  -5.5080732424212877 was WRONG. BUG: step_06_gold14.py/06b removed 1e-9 slivers at the VERTEX endpoints
  of the chamber [pi/2+1e-9, pi/2+pi/7-1e-9]; the p-ball integrand has a boundary layer there
  (|cos(th-pi/2)|^p amplifies p-derivatives by ~ln(1e-9)=-20.7), so the sliver error has err'''(1)~+1.05e-5.
  Being systematic (same for all p-samples), both Fourier radii + the FD stencil inherited the SAME err''',
  appearing self-consistent (7e-18) yet wrong. Verified by reproduction. 12-gon control passes.

## Task B — delta(N3): NONZERO (like octagon & 12-gon; not zero like decagon)
- Ported the 12-gon engine to F=Q(2cos(2pi/7)) exact arithmetic (engine14.py; FT arithmetic, Jloc/Jtab tables,
  pair/triple recursions). Validated vs quadrature to ~1e-59 on 13 integral types.
- mons14 (811 reduced monomials, weights {3:660,2:134,1:17}) -> assemble S (132 terms); **84*L*S = P'''(1)
  to 4.4e-50** (end-to-end validation of the whole engine+assembly).
- Weight-3 non-local symbol sigma(S): 8 terms, all in (span{alpha,u})^3 (x_i,1-x_i have no factor of 2).
  Atom symbols: L3x=li3, Ax=Li3(x)-Li3(1-x)+zeta(3)+prod (validated 1e-61), T3p=-(1/x)N repair (validated
  1e-61), T3a via S3 on verified G-words, Lam via DOD3's S_Li21 (-tau to PGEN1R convention).
- sigma(S) reduces to the 24-CR Li3 basis mod Decomp with **EXACT residual 0** (all 3 F-components):
  alphabet rule confirmed for the 14-gon. rank[li3|Decomp]=21, li3-span mod Decomp = 2-dim.
- delta(sigma) = sum c_z delta_li3(z): **delta vanishes on all 22 li3-relations** (well-defined); identical
  across 4+ particular solutions AND across a minimal 2-element basis AND both Lam conventions.
- RESULT: delta(N3) = c_a (a^u)(x)a + c_u (a^u)(x)u, with
    c_a = (601 - 182 a - 243 a^2)/32   ~ -0.118909811681569
    c_u = (613 - 118 a - 292 a^2)/32   ~ +0.369019755994578
  NONZERO. => 14-gon has genuine depth-2 (Li_{2,1}) trilogarithm content.
- C7 CHECK PASSES: the local zeta(3) part of S gives 84*L*z3_local = C_7 = -0.6838425376363484 to 3.6e-69,
  and the zeta(3) coefficient is ENTIRELY local (single/pair/triple families contribute 0).

## Task C — closed form of P'''(1): OPEN (not attempted; δ was the goal and is done)
"""
open('summary.md','w').write(summary)
print("summary.md updated")
