# s16: the edge law in terms of (edge length L, the two angles psi0, psi1 of the triangle O V0 V1 at V0 and V1,
# and the local weights w0, w1 of the two kink terms meeting the edge's endpoints, unit normals):
#   Z_e = 6 w0 w1 L^3 sin^3(psi0) sin^3(psi1) / sin^2(psi0+psi1) * (3 L sin(psi0) sin(psi1) (w0 cos psi0 + w1 cos psi1)/sin(psi0+psi1) - 2)
import mpmath as mp
from edge_num import edges, zeta3_edge_law
mp.mp.dps = 40
for name, bs, ws in (('hexA_w123', [mp.mpf(0), -mp.pi/3, -2*mp.pi/3], [1, 2, 3]),
                     ('hexB_angles', [mp.mpf(0), mp.mpf(1), mp.mpf('2.2')], [1, 1, 1]),
                     ('generic octagon', [mp.mpf(0), mp.mpf('0.7'), mp.mpf('1.6'), mp.mpf('2.5')], [1, mp.mpf('1.5'), mp.mpf('0.8'), mp.mpf('1.2')])):
    Z, parts = zeta3_edge_law(bs, ws); worst = 0
    for e, Ze in zip(edges(bs, ws, [0]*len(bs)), parts):
        V0, V1, L = e['V0'], e['V1'], e['L']
        d = (V1 - V0)/L
        r0 = mp.sqrt(V0[0]**2 + V0[1]**2); r1 = mp.sqrt(V1[0]**2 + V1[1]**2)
        psi0 = mp.acos(-(V0[0]*d[0] + V0[1]*d[1])/r0); psi1 = mp.acos((V1[0]*d[0] + V1[1]*d[1])/r1)
        w0, w1 = ws[e['j0']], ws[e['j1']]
        Zf = 6*w0*w1*L**3*mp.sin(psi0)**3*mp.sin(psi1)**3/mp.sin(psi0 + psi1)**2*(
             3*L*mp.sin(psi0)*mp.sin(psi1)*(w0*mp.cos(psi0) + w1*mp.cos(psi1))/mp.sin(psi0 + psi1) - 2)
        worst = max(worst, abs(Zf - Ze), abs(e['alpha'] - w0*L*mp.sin(psi0)), abs(e['beta'] - w1*L*mp.sin(psi1)))
    print(f"{name}: max |angle form - Z_e|, |alpha - w0 L sin psi0|, |beta - w1 L sin psi1| = {mp.nstr(worst, 3)}")
