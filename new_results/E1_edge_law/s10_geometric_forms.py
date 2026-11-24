# s10: geometric forms of the edge law, checked numerically.
#  Form 1: Z_e = 6 al be (h^2/L) (3(al u + be v) - 2), al = w_{j0} dist(V1, line OV0), be = w_{j1} dist(V0, line OV1),
#          u = |V0 F|/L, v = |F V1|/L (F = foot of the perpendicular from O), u + v = 1 (signed).
#  Form 2: al = L sin(eps0) / (2 h_prev), be = L sin(eps1) / (2 h_next)   (eps = exterior angle, h_prev/next =
#          distance from O to the neighbouring edge line at V0 / V1).  Pure geometry of K_1.
#  Form 3: tangential K_1 (all h_e equal): Z_e = (3/2) L sin(eps0) sin(eps1) (3 sin^2(eps0/2) + 3 sin^2(eps1/2) - 2).
import sys, mpmath as mp
sys.path.insert(0, 'research/08_theory/D15c_edge_zeta3')
import sympy as sp, sym_G
from edge_num import edges, zeta3_edge_law
mp.mp.dps = 40
G21, G12, S = sym_G.build()
args = [S[k] for k in ('s0', 's1', 'c0', 'c1', 'S0', 'Sc')]
fG = sp.lambdify(args, G21 + G12, 'mpmath')

def from_vertices(V):
    """centrally symmetric polygon with vertices V[0..m-1] (half, ccw, angles in [0,pi)) -> kinks b_j, weights w_j"""
    m = len(V)
    bs = [mp.atan2(v[1], v[0]) for v in V]
    A = mp.matrix(m, m)
    for k in range(m):
        for j in range(m):
            A[k, j] = abs(mp.sin(bs[j])*V[k][0]*(-1) + mp.cos(bs[j])*V[k][1])
    w = mp.lu_solve(A, mp.matrix([1]*m))
    return bs, [w[j] for j in range(m)]

def check(name, bs, ws, tangential=False):
    Es = edges(bs, ws, [0]*len(bs)); n = len(Es)
    Z, parts = zeta3_edge_law(bs, ws)
    w2 = mp.mpf(0); w3 = mp.mpf(0); wG = mp.mpf(0); wd = mp.mpf(0)
    for k, e in enumerate(Es):
        ep, en = Es[(k-1) % n], Es[(k+1) % n]
        def ext(e1, e2):  # exterior angle between direction of e1 and of e2
            d1 = (e1['V1'] - e1['V0'])/e1['L']; d2 = (e2['V1'] - e2['V0'])/e2['L']
            return mp.acos(d1[0]*d2[0] + d1[1]*d2[1])
        eps0, eps1 = ext(ep, e), ext(e, en)
        al2 = e['L']*mp.sin(eps0)/(2*ep['h']); be2 = e['L']*mp.sin(eps1)/(2*en['h'])
        w2 = max(w2, abs(al2 - e['alpha']), abs(be2 - e['beta']))
        r0 = mp.sqrt(e['V0'][0]**2 + e['V0'][1]**2)
        wd = max(wd, abs(e['alpha'] - ws[e['j0']]*e['L']*e['h']/r0))
        if tangential:
            Z3 = mp.mpf(3)/2*e['L']*mp.sin(eps0)*mp.sin(eps1)*(3*mp.sin(eps0/2)**2 + 3*mp.sin(eps1/2)**2 - 2)
            w3 = max(w3, abs(Z3 - parts[k]))
        a, b = e['th0'], e['th1']
        if b < a: b += 2*mp.pi
        def G(th):
            sv = [w*abs(mp.sin(th - bb)) for bb, w in zip(bs, ws)]
            cv = [w*mp.cos(th - bb)*mp.sign(mp.sin(th - bb)) for bb, w in zip(bs, ws)]
            return fG(sv[e['j0']], sv[e['j1']], cv[e['j0']], cv[e['j1']], mp.fsum(sv), mp.fsum(cv))
        wG = max(wG, abs(12*mp.quad(G, [a, b]) - parts[k]))
    print(f"{name}: Z = {mp.nstr(Z, 25)} | Form1 vs D15c-G integral {mp.nstr(wG, 3)} | al=w*L*h/|OV0| {mp.nstr(wd, 3)}"
          f" | Form2 al,be {mp.nstr(w2, 3)}" + (f" | Form3 {mp.nstr(w3, 3)}" if tangential else ''))
    return Es, parts

check('hexA_w123', [mp.mpf(0), -mp.pi/3, -2*mp.pi/3], [1, 2, 3])
check('hexB_angles', [mp.mpf(0), mp.mpf(1), mp.mpf('2.2')], [1, 1, 1])
check('regular_m7', [-(i*mp.pi/7) for i in range(7)], [1]*7, tangential=True)
# tangential non-regular hexagon: tangent lines to the unit circle with outer normals at angles 0.3, 1.4, 2.2 (+pi)
phis = [mp.mpf('0.3'), mp.mpf('1.4'), mp.mpf('2.2'), mp.mpf('0.3') + mp.pi]
V = []
for k in range(3):
    a1, a2 = phis[k], phis[k+1]
    M = mp.matrix([[mp.cos(a1), mp.sin(a1)], [mp.cos(a2), mp.sin(a2)]])
    v = mp.lu_solve(M, mp.matrix([1, 1])); V.append((v[0], v[1]))
bs, ws = from_vertices(V)
print("tangential hexagon: kinks", [mp.nstr(b, 8) for b in bs], "weights", [mp.nstr(w, 8) for w in ws])
Es, parts = check('tangential_hex', bs, ws, tangential=True)
print("   h_e:", [mp.nstr(e['h'], 12) for e in Es])
# tangential octagon (m = 4), normals at 0.1, 0.9, 1.5, 2.6
phis = [mp.mpf(x) for x in ('0.1', '0.9', '1.5', '2.6')] + [mp.mpf('0.1') + mp.pi]
V = []
for k in range(4):
    a1, a2 = phis[k], phis[k+1]
    M = mp.matrix([[mp.cos(a1), mp.sin(a1)], [mp.cos(a2), mp.sin(a2)]])
    v = mp.lu_solve(M, mp.matrix([1, 1])); V.append((v[0], v[1]))
bs, ws = from_vertices(V)
print("tangential octagon: weights", [mp.nstr(w, 8) for w in ws])
check('tangential_oct', bs, ws, tangential=True)
