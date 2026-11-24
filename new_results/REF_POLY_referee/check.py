import mpmath as mp, json, sys
from polys import polys
from geom import *
mp.mp.dps = 100
A = json.load(open(sys.argv[1])); B = json.load(open(sys.argv[2]))
z3, pi2, l2 = mp.zeta(3), mp.pi**2, mp.log(2)
def dig(x, y):
    e = abs(mp.mpf(x)-mp.mpf(y)); return '%.1e' % float(e) if e > 0 else '0'
for name, U in polys().items():
    a = [mp.mpf(x) for x in A[name]]; b = [mp.mpf(x) for x in B[name]]
    print('==', name)
    print('  run1 vs run2 |diff| for P,P1,P2,P3:', [dig(a[k], b[k]) for k in range(4)])
    P = b
    E = edges(U); Pv = mp.fsum(edge_data(U, *e)[0] for e in E)
    Ecal = mp.fsum(edge_data(U, *e)[1]**2/edge_data(U, *e)[0] for e in E)
    Zs = mp.fsum(Z_e(U, *e) for e in E)
    print('  P(1) - polygon perimeter:', dig(P[0], Pv), ' #edges', len(E))
    print("  P'(1) - referee general first-order formula:", dig(P[1], P1_general(U)))
    print("  P'(1) =", mp.nstr(P[1], 45))
    print('  E =', mp.nstr(Ecal, 30), ' sum Z_e =', mp.nstr(Zs, 30))
    if name.startswith('reg') or name == 'diamond':
        m = len(U); V = mp.cot(mp.pi/(2*m))
        R1 = mp.fsum(mp.sin(i*mp.pi/m)*mp.log(V/mp.sin(i*mp.pi/m)) for i in range(1, m))/V
        P1f = 4*m*(mp.sec(mp.pi/(2*m))-mp.cos(mp.pi/(2*m)))
        print('  thm:poly: P(1) formula diff', dig(P[0], P1f), "; P'(1)-P(1)R1:", dig(P[1], P1f*R1))
        s = mp.sin(mp.pi/(2*m)); c = mp.cos(mp.pi/(2*m))
        Cm = 48*m*s**4*c*(3*s**2-1)
        print('  sum Z_e - C_m:', dig(Zs, Cm), ' C_m =', mp.nstr(Cm, 30))
    if name in ('diamond', 'parallelogram'):
        print("  parallelogram law: P''(1)-pi^2E/3:", dig(P[2], pi2*Ecal/3), "; P'''(1)-6E(z3-pi^2/3):", dig(P[3], 6*Ecal*(z3-pi2/3)))
    if name in ('reg_hexagon', 'affreg_hexagon'):
        # 2^{d/(1+d)} [P0 + pi^2 E/24 d^2 - E/8 (z3 + 2pi^2/3) d^3]
        g = [0, l2, -l2+l2**2/2, l2-l2**2+l2**3/6]   # exp(l2*(d-d^2+d^3)) coefficients
        f = [Pv, 0, pi2*Ecal/24, -Ecal/8*(z3+2*pi2/3)]
        ser = [mp.fsum((1 if k == 0 else g[k])*0 for k in range(1)) for _ in range(4)]
        gg = [mp.mpf(1)] + g[1:]
        ser = [mp.fsum(gg[k]*f[n-k] for k in range(n+1)) for n in range(4)]
        D = [ser[0], ser[1], 2*ser[2], 6*ser[3]]
        print('  hexagon law diffs P..P3:', [dig(P[k], D[k]) for k in range(4)])
        print('  coefficient of zeta(3) in P\'\'\'(1) from law: -3E/4 =', mp.nstr(-3*Ecal/4, 25), ' sum Z_e =', mp.nstr(Zs, 25))
    if name == 'reg_hexagon':
        pr = mp.sqrt(3)/8*(16*l2**3+3*pi2*l2-96*l2**2+96*l2-9*z3-6*pi2)
        print("  printed regular hexagon P'''(1) diff:", dig(P[3], pr), "; P'(1)-2sqrt3 log2:", dig(P[1], 2*mp.sqrt(3)*l2))
