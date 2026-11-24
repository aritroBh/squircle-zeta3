# g9_direct5.py -- fresh numerical check of P^(5)(1), by finite differences in p of the perimeter of K_p (own code).
# K_p = {sum_j w_j^{1-p} |<u_j,x>|^p <= 1}; K_1 = {sum_j |<u_j,x>| <= 1}.  The perimeter is the length of the curve
# x_p(t) = y(t) F_p(y(t))^{-1/p} over every edge y(t) = V0 + t (V1 - V0) of K_1 (Lemma g4:lem:arc); no series in p is used.
# P^(5)(1) = 5-point-pair central difference (steps h, 2h, 3h) with Richardson extrapolation in h^2.
import mpmath as mp, sys, itertools
def corners_from_u(us):
    m = len(us); pts = []
    for (ux, uy) in us:
        dx, dy = -uy, ux
        nrm = sum(abs(ux*dx_ + uy*dy_) for (ux, uy) in us for (dx_, dy_) in [(dx, dy)])
        for sg in (1, -1): pts.append((sg*dx/nrm, sg*dy/nrm))
    pts.sort(key=lambda P: mp.atan2(P[1], P[0]))
    return pts
def us_from_corners(C, mode='half'):
    """C: corners C_0..C_{m-1} (consecutive, counter-clockwise, half of the polygon). Returns u_j (kink line through C_j), via Lemma g4:lem:polygon(b)."""
    m = len(C); allC = list(C) + [(-x, -y) for (x, y) in C]
    def grad(Pa, Pb):               # g with <g,Pa> = <g,Pb> = 1
        det = Pa[0]*Pb[1] - Pa[1]*Pb[0]
        return ((Pb[1]-Pa[1])/det, (Pa[0]-Pb[0])/det)
    g = [grad(allC[k], allC[(k+1) % (2*m)]) for k in range(2*m)]
    us = []
    for k in range(m):
        gk, gp = g[k], g[k-1]
        us.append(((gk[0]-gp[0])/2, (gk[1]-gp[1])/2))
    return us
def perimeter(p, us, ws, edges):
    d = p - 1
    tot = mp.mpf(0)
    for (V0, V1) in edges:
        yp = (V1[0]-V0[0], V1[1]-V0[1])
        def f(t):
            y = (V0[0]+t*yp[0], V0[1]+t*yp[1])
            F = mp.mpf(0); D = mp.mpf(0)
            for (u, w) in zip(us, ws):
                s = u[0]*y[0] + u[1]*y[1]; a = abs(s)
                ap = (u[0]*yp[0] + u[1]*yp[1])*(1 if s >= 0 else -1)
                fac = (a/w)**d
                F += a*fac; D += ap*fac
            Phi = D/F
            vx, vy = yp[0] - Phi*y[0], yp[1] - Phi*y[1]
            return F**(-1/p)*mp.sqrt(vx*vx + vy*vy)
        tot += mp.quad(f, [0, mp.mpf(1)/2, 1])
    return tot
def P5(us, ws, h0=mp.mpf('0.01'), levels=5):
    C = corners_from_u(us)
    m = len(us)
    edges = [(C[k], C[(k+1) % (2*m)]) for k in range(2*m)]
    # use central symmetry: edges k and k+m have equal length elements
    half = edges[:m]
    def P(p): return 2*perimeter(p, us, ws, half)
    rows = []
    for l in range(levels):
        h = h0/2**l
        v = (P(1+3*h) - 4*P(1+2*h) + 5*P(1+h) - 5*P(1-h) + 4*P(1-2*h) - P(1-3*h))/(2*h**5)
        rows.append(v)
    T = [rows]
    for k in range(1, levels):
        prev = T[-1]
        T.append([(4**k*prev[i+1] - prev[i])/(4**k - 1) for i in range(len(prev)-1)])
    return T, P(mp.mpf(1))
def hexagon(Vs, ws):
    us = us_from_corners(Vs)
    return us, ws
if __name__ == '__main__':
    mp.mp.dps = int(sys.argv[4]) if len(sys.argv) > 4 else 60
    which = sys.argv[1]
    if which == 'reg6':
        # regular hexagon, K_p = sum |<n_i,x>|^p <= 1
        us = [(mp.cos(i*mp.pi/3), mp.sin(i*mp.pi/3)) for i in (1, 2, 3)]; ws = [1, 1, 1]
        known = mp.mpf('656.41536196762475447664760905831120964826426464898548604185355096214310290964572236224057164257468405023297944386495')
    elif which == 'golden':
        phi = (1+mp.sqrt(5))/2
        us = us_from_corners([(mp.mpf(1), mp.mpf(0)), (phi, mp.mpf(1)), (mp.mpf(0), mp.mpf(1))]); ws = [2, 3, 5]
        known = mp.mpf('-78.358599873964602127531124178168617182820')
    elif which == 'paper':
        us = us_from_corners([(mp.mpf(1), mp.mpf(0)), (mp.mpf(2), mp.mpf(2)), (mp.mpf(0), mp.mpf(1))]); ws = [1, 3, 1]
        known = mp.mpf('-449.23396800163978089849025140306005176022')
    elif which == 'shear':
        us = us_from_corners([(mp.mpf(2), mp.mpf(0)), (mp.mpf(1), mp.mpf(2))]); ws = [1, 1]
        z2, z3, z5 = mp.zeta(2), mp.zeta(3), mp.zeta(5); s5, s13 = mp.sqrt(5), mp.sqrt(13)
        known = (z5*(-6912*s5/25 - 172800*s13/28561) + z2*z3*(-10752*s5/25 - 821760*s13/28561) + z2**2*(119808*s5/125 + 1262592*s13/28561)
                 + z3*(4608*s5/5 + 23040*s13/169) + z2*(-3072*s5/5 - 15360*s13/169))
    elif which == 'oct':
        us = [(mp.cos(i*mp.pi/4), mp.sin(i*mp.pi/4)) for i in (1, 2, 3, 4)]; ws = [1, 1, 1, 1]
        known = mp.mpf('164.12775911330597065104567331245436130507711746936')
    h0 = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('0.01'); levels = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    print('dps', mp.mp.dps, 'h0', h0, 'levels', levels)
    T, P1 = P5(us, ws, h0, levels)
    print(which, ' P(1) =', mp.nstr(P1, 20))
    for k, row in enumerate(T):
        print('  Richardson level', k, ':', [mp.nstr(v, 25) for v in row][-2:])
    best = T[-1][-1]
    print('  P5(1) fd  =', mp.nstr(best, 30))
    print('  reference =', mp.nstr(known, 30))
    print('  |difference| =', mp.nstr(abs(best-known), 3))
