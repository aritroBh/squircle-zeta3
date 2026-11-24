"""Exact (instance-wise) derivation of the parallelogram and affine-regular-hexagon laws from the
referee's exact edge expansion: rational corners, exact Beta table, exact comparison."""
import sympy as sp, pickle, itertools, sys
from edge_sym import t, L, al, be, B1, C1, x0, y0, x1, y1, X0, X1, Xf1, Xf2
from exact_laws_lib import mom, lt, l1, Bfun, xx, yy
E = pickle.load(open('edge_E_general.pkl', 'rb'))
print('Beta lemma check:'); ok = True
for n, k in itertools.product(range(4), range(3)):
    for a, b in [(3,0),(0,3),(2,1),(1,2),(2,0),(1,1),(0,2),(1,0),(0,1)]:
        v = sp.expand(sp.simplify(sp.expand_func(sp.diff(Bfun, xx, a, yy, b).subs({xx: n+1, yy: k+1}))))
        c3 = v.coeff(sp.zeta(3))
        want = 2*sp.gamma(n+1)*sp.gamma(k+1)/sp.gamma(n+k+2) if (a, b) in [(2,1),(1,2)] else 0
        if sp.simplify(c3-want) != 0 or sp.expand(v-c3*sp.zeta(3)).has(sp.zeta(3)): ok = False; print(' FAIL', n, k, a, b)
print('  108 cases as claimed:', ok, flush=True)
l2 = sp.log(2)
def edge_int(k, V0, V1, subs):
    Lv = sp.sqrt((V1[0]-V0[0])**2+(V1[1]-V0[1])**2)
    e = sp.expand(E[k].subs(subs).subs({x0: V0[0], y0: V0[1], x1: V1[0], y1: V1[1], L: Lv}))
    P = sp.Poly(e, t, lt, l1)
    return sp.expand(sum(c*mom(n, a, b) for (n, a, b), c in P.terms()))
def hL(V0, V1):
    Lv = sp.sqrt((V1[0]-V0[0])**2+(V1[1]-V0[1])**2)
    cr = V0[0]*V1[1]-V0[1]*V1[0]
    return cr**2/Lv**3, Lv
Q = sp.Rational
sub_par = {al: 1, be: 1, B1: 0, C1: 0, X0: lt, X1: l1, Xf1: 0, Xf2: 0}
for V0, V1 in [((Q(1),Q(0)),(Q(1,3),Q(2))), ((Q(2),Q(-1,2)),(Q(-1,5),Q(3,2)))]:
    eds = [(V0, V1), (V1, (-V0[0], -V0[1]))]
    Ecal = 2*sum(hL(*e)[0] for e in eds); P1 = 2*sum(hL(*e)[1] for e in eds)
    law = [P1, 0, sp.pi**2*Ecal/6, Ecal*(sp.zeta(3)-sp.pi**2/3)]
    for k in range(4):
        tot = 2*sum(edge_int(k, *e, sub_par) for e in eds)
        print('parallelogram', V0, V1, 'd^%d: total - law =' % k, sp.simplify(sp.radsimp(sp.expand(tot-law[k]))), flush=True)
sub_hex = {al: Q(1,2), be: Q(1,2), B1: Q(1,2), C1: 0, X0: lt-l2, X1: l1-l2, Xf1: -l2, Xf2: 0}
g = [1, l2, -l2+l2**2/2, l2-l2**2+l2**3/6]
for V0, V1 in [((Q(1,2),Q(0)),(Q(1,2),Q(1,2))), ((Q(1),Q(1,3)),(Q(2,5),Q(3,2)))]:
    V2 = (V1[0]-V0[0], V1[1]-V0[1])
    eds = [(V0, V1), (V1, V2), (V2, (-V0[0], -V0[1]))]
    Ecal = 2*sum(hL(*e)[0] for e in eds); P1 = 2*sum(hL(*e)[1] for e in eds)
    f = [P1, 0, sp.pi**2*Ecal/24, -Ecal/8*(sp.zeta(3)+2*sp.pi**2/3)]
    for k in range(4):
        law = sum(g[j]*f[k-j] for j in range(k+1))
        tot = 2*sum(edge_int(k, *e, sub_hex) for e in eds)
        print('aff-reg hexagon', V0, V1, 'd^%d: total - law =' % k, sp.simplify(sp.radsimp(sp.expand(tot-law))), flush=True)
