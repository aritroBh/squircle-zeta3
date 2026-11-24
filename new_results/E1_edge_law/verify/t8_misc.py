# T8: (1) FD-only total of the rational hexagon vs its PSLQ closed form; (2) affine-regular hexagon law on the verifier's own
# non-regular affine-regular hexagon V=(1,0),(1,2),(0,2) (power convention), by the jet code at 60 digits;
# (3) C_m asymptotics, numerically and exactly.
import sys, json, mpmath as mp, sympy as sp
sys.path.insert(0, 'research/10_frontier/E1_edge_law/verify')
from vjet import *
mp.mp.dps = 60
z3, l2, pi2 = mp.zeta(3), mp.log(2), mp.pi**2
s2, s5 = mp.sqrt(2), mp.sqrt(5)
closed = s5*(mp.mpf(231)/125*z3 + mp.mpf(111)/125*pi2*l2 + mp.mpf(52)/125*l2**3 - mp.mpf(16)/25*pi2 - mp.mpf(288)/25*l2**2 + mp.mpf(96)/5*l2) \
       + s2*(-mp.mpf(15)/64*z3 + pi2*l2/16 + 16*l2**3 - pi2/16 - 48*l2**2 + 24*l2)
fd = json.load(open('research/10_frontier/E1_edge_law/verify/out/t4_fd_170_0.0004_26.json'))
print("(1) rational hexagon: FD-only total P'''(1) =", mp.nstr(mp.mpf(fd['total'][3]), 50))
print("    PSLQ closed form                      =", mp.nstr(closed, 50), " |diff| =", mp.nstr(abs(mp.mpf(fd['total'][3]) - closed), 3))
print("    closed form with the EDGE-LAW zeta(3) coefficient instead differs by (159 sqrt5/125) z3 =", mp.nstr(159*s5/125*z3, 20))
# (2) affine-regular hexagon, not regular, not H_p
V = [(1, 0), (1, 2), (0, 2)]
vs = [(-mp.mpf(v[1])/4, mp.mpf(v[0])/4) for v in V]     # <v_j,x> = cross(V_j,x)/4
cs = [mp.mpf(1)]*3; lams = [0]*3
ks = kink_angles(vs); pts = ks + [ks[0] + 2*mp.pi]
tot = [mp.mpf(0)]*4
for a, b in zip(pts[:-1], pts[1:]):
    c = arc_coeffs_all(a, b, vs, cs, lams); tot = [x + y for x, y in zip(tot, c)]
P = [tot[0], tot[1], 2*tot[2], 6*tot[3]]
S = 9 + 8/(5*s5); P1 = 6 + 2*s5; lam = l2
law = [P1, P1*lam, 2*(P1*(lam**2/2 - lam) + pi2*S/24), 6*(P1*(lam**3/6 - lam**2 + lam) + lam*pi2*S/24 - S/8*(z3 + 2*pi2/3))]
print("(2) affine-regular hexagon V=(1,0),(1,2),(0,2), power convention, S = 9 + 8/(5 sqrt5):")
for k in range(4):
    print(f"    P^({k})(1): jet {mp.nstr(P[k], 30)}  law {mp.nstr(law[k], 30)}  |diff| {mp.nstr(abs(P[k]-law[k]), 3)}")
print("    zeta(3) coefficient of the law = -(3/4) S =", mp.nstr(-3*S/4, 25))
# (3) C_m
m = sp.symbols('m', positive=True)
x = sp.pi/(2*m)
Cm = 48*m*sp.sin(x)**4*sp.cos(x)*(3*sp.sin(x)**2 - 1)
ser = sp.series(Cm.subs(m, 1/sp.Symbol('e')), sp.Symbol('e'), 0, 9).removeO()
print("(3) C_m in powers of e = 1/m:", sp.expand(ser))
for mm in (10, 100, 1000):
    v = sp.N(Cm.subs(m, mm), 20)
    print(f"    m={mm}: m^3 C_m = {sp.N(v*mm**3, 15)}; m^3(-3pi^4/m^3 + 25pi^6/(8m^5)) = {sp.N(-3*sp.pi**4 + 25*sp.pi**6/(8*mm**2), 15)}; old ledger 9pi^5/(16 m^4)*m^3 = {sp.N(9*sp.pi**5/(16*mm), 6)}")
