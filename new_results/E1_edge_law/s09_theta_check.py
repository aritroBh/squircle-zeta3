# s09: link to D15c's polar-angle (theta) framework.
# (i) pointwise D_k(theta) dtheta/dt = E_k(t) (same curve, p-independent reparametrisation theta(t));
# (ii) D15c's symbolic G21, G12 (sym_G.py), made weighted (s_i -> w_i s_i, c_i -> w_i c_i), satisfy
#      G21(theta(t)) dtheta/dt = g21(t), G12 dtheta/dt = g12(t) pointwise (uniqueness of the log^2 t log(1-t)
#      coefficient), hence int_arc (G21+G12) dtheta = int_0^1 (g21+g12) dt = Z_e/12.
import sys, pickle, mpmath as mp, sympy as sp
sys.path.insert(0, 'research/08_theory/D15b_polygon_law')
sys.path.insert(0, 'research/08_theory/D15c_edge_zeta3')
import poly_coeffs as pc, sym_G
from edge_num import edges, E_series, zeta3_edge_law
mp.mp.dps = 40
G21, G12, S = sym_G.build()
args = [S[k] for k in ('s0', 's1', 'c0', 'c1', 'S0', 'Sc')]
print("free symbols of D15c G21:", sorted(map(str, G21.free_symbols)))
fG21 = sp.lambdify(args, G21, 'mpmath'); fG12 = sp.lambdify(args, G12, 'mpmath')
D = pickle.load(open('s01_E.pkl', 'rb'))
t, al, be, L, h, s0 = sp.symbols('t alpha beta L h s0', real=True)
fg21 = sp.lambdify((t, al, be, L, h, s0), D['g21'], 'mpmath'); fg12 = sp.lambdify((t, al, be, L, h, s0), D['g12'], 'mpmath')
for name, bs, ws, lams in (('hexA_w123', [mp.mpf(0), -mp.pi/3, -2*mp.pi/3], [1, 2, 3], [0, 0, 0]),
                           ('hexB_angles', [mp.mpf(0), mp.mpf(1), mp.mpf('2.2')], [1, 1, 1], [0, 0, 0]),
                           ('reg_m5', [-(i*mp.pi/5) for i in range(5)], [1]*5, [0]*5)):
    coeffs = pc.make_coeffs(bs, ws, lams)
    Z, parts = zeta3_edge_law(bs, ws)
    wD = mp.mpf(0); wG = mp.mpf(0); wI = mp.mpf(0)
    for e, Ze in zip(edges(bs, ws, lams), parts):
        V0, V1 = e['V0'], e['V1']
        def theta_of(tt):
            y = (1-tt)*V0 + tt*V1
            th = mp.atan2(y[1], y[0])
            if th < e['th0'] - mp.mpf('1e-20'): th += 2*mp.pi
            return th, y
        def wG_at(th, y):
            sv = [w*abs(mp.sin(th - b)) for b, w in zip(bs, ws)]
            cv = [w*mp.cos(th - b)*mp.sign(mp.sin(th - b)) for b, w in zip(bs, ws)]
            a_ = [sv[e['j0']], sv[e['j1']], cv[e['j0']], cv[e['j1']], mp.fsum(sv), mp.fsum(cv)]
            return fG21(*a_), fG12(*a_)
        for tt in (mp.mpf('0.13'), mp.mpf('0.5'), mp.mpf('0.77')):
            th, y = theta_of(tt)
            jac = e['h']*e['L']/(y[0]**2 + y[1]**2)
            Dk = coeffs(th); Ek = E_series(e, tt)
            wD = max(wD, max(abs(Dk[k]*jac - Ek[k]) for k in range(4)))
            g21t = fg21(tt, e['alpha'], e['beta'], e['L'], e['h'], e['s0'])
            g12t = fg12(tt, e['alpha'], e['beta'], e['L'], e['h'], e['s0'])
            G21v, G12v = wG_at(th, y)
            wG = max(wG, abs(G21v*jac - g21t), abs(G12v*jac - g12t))
        a, b = e['th0'], e['th1']
        if b < a: b += 2*mp.pi
        IG = mp.quad(lambda th: sum(wG_at(th, None)), [a, b])
        wI = max(wI, abs(12*IG - Ze))
    print(f"{name}: max pointwise |D_k dth/dt - E_k| = {mp.nstr(wD, 3)};  max |G dth/dt - g| = {mp.nstr(wG, 3)};"
          f"  max per edge |12 int(G21+G12)dth - Z_e| = {mp.nstr(wI, 3)};  Z = {mp.nstr(Z, 25)}")
