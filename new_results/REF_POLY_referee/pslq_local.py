"""Integer relations for the fully local polygons: is the zeta(3) coefficient of P'''(1) equal to sum Z_e?"""
import mpmath as mp
from jet import perimeter_jet
from geom import edges, Z_e, edge_data
mp.mp.dps = 140
R = lambda a, b=1: mp.mpf(a)/b
cases = {
  'parallelogram (1,1/2),(-1/3,1)': [(R(1),R(1,2)),(R(-1,3),R(1))],
  'regular hexagon': [(mp.cos(i*mp.pi/3), mp.sin(i*mp.pi/3)) for i in (1,2,3)],
  'affine-regular hexagon |x|+|y|+|x-y|': [(R(1),R(0)),(R(0),R(1)),(R(1),R(-1))],
}
z3, pi2, l2 = mp.zeta(3), mp.pi**2, mp.log(2)
import sys
for name, U in cases.items():
    if len(sys.argv)>1 and not name.startswith(sys.argv[1]): continue
    P = perimeter_jet(U, h=mp.mpf(1)/256)
    P3 = 6*P[3]
    Zs = mp.fsum(Z_e(U, *e) for e in edges(U))
    Ls = sorted(set(mp.nstr(edge_data(U, *e)[0]**2, 40) for e in edges(U)))
    print('==', name, " P'''(1) =", mp.nstr(P3, 60)); print('   squared edge lengths:', Ls, ' sum Z_e =', mp.nstr(Zs, 30))
    if name.startswith('parallelogram'):
        # L^2 = 5/... ; surds
        Lsq = sorted(set([edge_data(U, *e)[0]**2 for e in edges(U)]), key=float)
        dd = []
        for x in Lsq:
            if all(abs(x-y) > 1e-50 for y in dd): dd.append(x)
        Lsq = dd
        s1, s2 = [mp.sqrt(x) for x in Lsq]
        print('   L^2 identified:', [mp.identify(x) for x in Lsq])
        basis = [P3] + [c*s for s in (1/s1, 1/s2) for c in (z3, pi2, 1)]
        names = ['P3'] + ['%s/L%d' % (c, i) for i in (1,2) for c in ('z3','pi2','1')]
    elif name == 'regular hexagon':
        basis = [P3/mp.sqrt(3), z3, pi2*l2, l2**3, pi2, l2**2, l2, 1]
        names = ['P3/sqrt3', 'z3', 'pi2*l2', 'l2^3', 'pi2', 'l2^2', 'l2', '1']
    else:
        basis = [P3] + [s*c for s in (1, mp.sqrt(2)) for c in (z3, pi2*l2, l2**3, pi2, l2**2, l2, 1)]
        names = ['P3'] + ['%s%s' % (s, c) for s in ('', 'r2*') for c in ('z3','pi2*l2','l2^3','pi2','l2^2','l2','1')]
    for dps in (90, 130):
        with mp.workdps(dps):
            rel = mp.pslq([+x for x in basis], maxcoeff=10**8, maxsteps=10**6)
        if rel is None: print('   dps', dps, ': no relation'); continue
        res = mp.fsum(r*x for r, x in zip(rel, basis))
        print('   dps', dps, 'relation:', list(zip(names, rel)), ' residual %.1e' % float(abs(res)))
        # zeta(3) coefficient
        if name.startswith('parallelogram'):
            c = -(rel[1]/s1 + rel[4]/s2)/rel[0]
        elif name == 'regular hexagon':
            c = -mp.sqrt(3)*rel[1]/rel[0]
        else:
            c = -(rel[1] + rel[8]*mp.sqrt(2))/rel[0]
        print('   zeta(3) coefficient from relation =', mp.nstr(c, 30), '; minus sum Z_e = %.1e' % float(abs(c-Zs)))
