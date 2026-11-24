# s15: independent numerical check (D15b polar coefficient integrals) of the two family laws, power convention:
#  parallelogram  sum_{j=1,2} |<v_j,x>|^p <= 1:     P = P1 + (pi^2 S/6) d^2 + S (zeta3 - pi^2/3) d^3 + O(d^4)
#  affine-regular hexagon T(regular), v_j = T^{-T} n_j:  P = 2^{1-1/p} [P1 + (pi^2 S/24) d^2 - (S/8)(zeta3 + 2pi^2/3) d^3] + O(d^4)
# with S = sum over edges of h_e^2/L_e.
import sys, mpmath as mp
sys.path.insert(0, 'research/08_theory/D15b_polygon_law')
import poly_coeffs as pc
from edge_num import edges
mp.mp.dps = 45
def conv(vs):
    bs, ws, lams = [], [], []
    for v in vs:
        nv = mp.sqrt(v[0]**2 + v[1]**2); n = (v[0]/nv, v[1]/nv)
        bs.append(mp.atan2(-n[0], n[1])); ws.append(nv); lams.append(mp.log(nv))
    return bs, ws, lams
def check(name, vs, hexa):
    bs, ws, lams = conv(vs)
    Es = edges(bs, ws, lams)
    S = mp.fsum([e['h']**2/e['L'] for e in Es]); P1 = mp.fsum([e['L'] for e in Es])
    c = pc.P_coeffs_general(bs, ws, lams, kmax=3)
    P = [c[0], c[1], 2*c[2], 6*c[3]]
    if not hexa:
        pred = [P1, 0, mp.pi**2*S/3, 6*S*(mp.zeta(3) - mp.pi**2/3)]
    else:
        a = mp.log(2); Q2 = mp.pi**2*S/12; Q3 = -mp.mpf(3)/4*S*(mp.zeta(3) + 2*mp.pi**2/3)
        pred = [P1, P1*a, Q2 + P1*(a**2 - 2*a), Q3 + 3*a*Q2 + P1*(6*a - 6*a**2 + a**3)]
    print(name, " S =", mp.nstr(S, 20), " alpha,beta on edges:", sorted(set(mp.nstr(e['alpha'], 10) for e in Es)),
          " |P^(k) - law| k=0..3:", [mp.nstr(abs(x - y), 3) for x, y in zip(P, pred)])
check("sheared parallelogram", [(mp.mpf(1), mp.mpf('-0.3')), (mp.mpf('0.5'), mp.mpf('1.2'))], False)
check("parallelogram 2", [(mp.mpf('0.7'), mp.mpf('0.2')), (mp.mpf('-0.4'), mp.mpf('1.9'))], False)
T = mp.matrix([[1, mp.mpf('0.4')], [mp.mpf('0.2'), mp.mpf('1.3')]]); Tit = (T**-1).T
vs = []
for j in range(3):
    n = mp.matrix([mp.cos(j*mp.pi/3), mp.sin(j*mp.pi/3)]); v = Tit*n; vs.append((v[0], v[1]))
check("affine-regular hexagon T(regular)", vs, True)
