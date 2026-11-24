import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
summary = r'''# AV2 summary (Beukers' zeta(3) integrals on squircle domains)

## Objects
- D_p = {x,y>=0: x^p+y^p<=1}, B_p octant of 3D superball.
- B2(p) = -(1/2)∬_{D_p} log(xy)/(1-xy)  (= zeta(3) at p=inf)
- B3(p) = ∭_{B_p} 1/(1-(1-xy)z)          (= 2 zeta(3) at p=inf)

## Proved reductions (in av2core.py)
- B2 = sum_{n>=1} phi(n/p)/n^3,  phi(w) = w^3 Gamma(w)^2/Gamma(2w+1) (psi(2w+1)-psi(w)), phi(0)=1.
- B2 = (1/p)∫_0^{pi/4} (cosθ sinθ)^{2/p-1} g(c) dθ, c=(cosθ sinθ)^{2/p}, g(c)=(Li_2(c)+ln c ln(1-c))/c.
  Layer-free v-form: (1/2)∫_0^∞ u^{2/p} g(u^{2/p}) dv, u=1/(2 cosh(pv/2)).
  c-form: (1/2)∫_0^{2^{-2/p}} g(c)/sqrt(1-4 c^p) dc.
- B3 = ∬_{D_p} -log(1-(1-xy)Z)/(1-xy), Z=(1-x^p-y^p)^{1/p}
     = 2∫_0^{c*} F(c,p)/sqrt(1-4c^p) dc, F(c,p)=∫_0^1 r(-log(1-(1-r^2 c)(1-r^p)^{1/p}))/(1-r^2 c) dr
     = peak-subtracted: 2F(c*) c* B(1/p,1/2)/p + 2∫[F(c)-F(c*)]/sqrt(1-4c^p) dc.

## Main results
- B2(p) table: 24 values, 110 digits, series+quad agree ~1e-109.
- B3(p) table: 19 values p=1..2^20, 45-130 digits (nested tanh-sinh, two outer levels + two inner grids).
- Square end (PROVED via Euler-Maclaurin):
  B2(p) = zeta(3) + C2 e^2 + (zeta3/2) e^3 - (3 zeta4/16) e^4 + 0 e^5 + (phi6/120) e^6 + 0 e^7 + ..., e=1/p,
  phi(w)=Taylor coeffs exact: phi_2=0, phi_3=-z3, phi_4=pi^4/40, phi_6=79 pi^6/7560 - 4 z3^2, ...
  C2 = ∫(phi-1)w^{-3}dw = -zeta(2)/2 = -pi^2/12  [NUMERICAL, 81 digits + 54-digit fit].
  NO log terms anywhere in B2.
- B3 square end: B3 = 2zeta3 + 2 zeta(2) ln(eps) eps^2 + A2 eps^2 + 4 zeta(3) eps^3 + ...
  ln-coefficient = 2 zeta(2) [NUMERICAL, 11 digits, mechanism = crossover at c~e^{-p}];
  A2 = -5.4189939975 [OPEN]; eps^3 coeff = 4 zeta(3) [NUMERICAL ~5 digits]. E = B3-2B2 studied separately.
- Triangle end: B2(1) = 2 zeta(3)/3 [NUMERICAL 101 digits + PSLQ; proof sketch].
  B3(1) = 0.246617302384923740565168068839592055997186101518 (123 digits) [OPEN; P1+P2 split].
  B2 tau-expansion: d1 = 0.524664448852056928181428100965551051111928417616155814312367 (exact derivative) [OPEN].
  B3 tau-expansion: d1 = 0.7266961 [NUMERICAL].
- B2(2): even-index part = zeta(3)/12 [NUMERICAL]; odd part = 0.95919692720272038 [OPEN];
  constituent sums S2 = Cl_2(pi/3)-1, S3 = (pi-3)/3 [exact, NUMERICAL 106 digits].
- Legendre (n<=12): p=inf Beukers verified: d_n^3 L_n = A_n zeta3 + B_n, A_n = 2 apery(n) d_n^3, B_n integer;
  |L_n|^{1/n} -> 0.0294 = (sqrt2-1)^4. Finite p: structure does NOT persist (r_n = L_n - 2 apery B2(p)
  not rational at any finite p tested); decay rates |L_n|^{1/n} grow from 0.0294 (inf) to ~0.47 (p=1, n=12);
  diagonal moments M(j,j,p) are distinct p-periods; at p=inf they collapse to 2(zeta3 - H_j^{(3)}).

## Files
- av2core.py: all formulas/methods. b2_table.json, b3_final.json, b3_hifi.json, b2_largep.json, c2.json,
  legendre_inf.json, legendre_finite.json.
- step_02..step_47 + b3run.py, b3run2.py: computation scripts (see progress.log).
'''
open('summary.md','w').write(summary)
with open('progress.log','a') as f:
    f.write("summary.md written; final report assembly. DONE\n")
print("summary.md written", len(summary))
