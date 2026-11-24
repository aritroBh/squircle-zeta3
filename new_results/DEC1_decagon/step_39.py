import os
os.chdir('research/08_theory/DEC1_decagon')
res = r"""# DEC1 — the regular decagon's P'''(1) in closed form (work in progress)

## Bottom line so far
- The decagon edge has THREE non-local terms, but by the reflection symmetry of the edge their ratios are
  {phi, 1, 1/phi}: **the middle non-local term is constant** (kappa=0). So genuinely there are TWO non-local
  logs, Landen-paired (x1 = -1/phi <-> x4 = 1/phi^2, x4 = x1/(x1-1)) — the same integral-type inventory as
  the octagon, plus one extra constant letter log(2 phi).
- Edge data (PROVED, exact): A0 t, A1(1-t) with A0=A1=1/(2 phi^2); B1(1+kap1 t), B1=1/(2 phi^2), kap1=phi-1=1/phi;
  B4(1+kap4 t), B4=1/(2 phi), kap4=1/phi-1=-1/phi^2; B5=1/(2 phi) constant. Field: Q(sqrt5) for all
  coefficients; the ONLY radical is the overall edge length L = sqrt((25-11 sqrt5)/10) (analog of the octagon's u).
- E3 built exactly: 55 log-monomials (5 logs L0..L4, weight <=3), coefficients (p+q sqrt5) t^n, n<=3, max
  denominator 384, reconstruction residual 1.6e-99. In the reduced basis {log t, log(1-t), w1=log(1+t/phi),
  w4=log(1-t/phi^2), l2=log 2, lf=log phi}: 261 monomials = 95 local + 70+70 single-ratio + 26 mixed
  (20 distinct mixed integrals: types (a,b,c1,c4) in {(0,0,1,1),(1,0,1,1),(0,1,1,1),(0,0,2,1),(0,0,1,2)},
  n = 0..3 — EXACTLY the octagon's five mixed types).
- P'''(1) = 3!·(10 edges)·int E3 = 60·L·int E3n.
- Local zeta(3) coefficient from E1's edge law: C_5 = 10 Z_e = sqrt(10-2 sqrt5)(105 sqrt5-255)/32
  = -1.48510280041350954453523545754949332249477… matches D15b's 120-digit PSLQ value (all 120 digits).

## Verification so far
| check | residual |
|---|---|
| E3(0.3) exact vs direct d-series (80 dps) | 6.3e-82 |
| 60 L int E3n vs 200-digit P'''(1) (60 dps quad) | 9.3e-59 |
| Q(sqrt5) reconstruction of coefficients | 1.6e-99, max denom 384 |

## To do
600-digit tanh-sinh P'''(1); Jloc/Jtab checks; mixed engine (D-constants); A-atom and depth-2 reductions
(Lam = Li_{2,1}(-1/phi,-1/phi) = Li_3(1/phi^2) by double-shuffle; swap = Li_{2,1}(1/phi^2,-phi)); assembly; conjecture.
"""
open('result.md','w').write(res)
with open('progress.log','a') as f: f.write("step: early result.md written.\n")
print("result.md written")