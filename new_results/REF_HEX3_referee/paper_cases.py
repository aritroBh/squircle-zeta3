"""Referee re-computation of the three test hexagons quoted after thm:hexagon, and of the affine-regular
hexagon law of prop:edgelaw (to which thm:hexagon is said to reduce at r=1)."""
import sympy as sp, mpmath as mp, direct, formula
from run_tests import frame
mp.mp.dps = 60
phi = (1+sp.sqrt(5))/2
cases = {
 'paper 1: (1,0),(3,3),(0,1) w=3,5,3': ([(1,0),(3,3),(0,1)], [3,5,3]),
 'paper 2: (2,0),(3,2),(0,1) w=2,5,7': ([(2,0),(3,2),(0,1)], [2,5,7]),
 'paper 3: golden (1,0),(phi,1),(0,1) w=2,3,5': ([(1,0),(phi,1),(0,1)], [2,3,5]),
 'affine-regular (frame a=b=1) M=[[2,1],[0,1]] w=1,1,1': (frame(1,1,[[2,1],[0,1]]), [1,1,1]),
 'affine-regular (frame a=b=1) M=[[2,1],[0,1]] w=2,3,5': (frame(1,1,[[2,1],[0,1]]), [2,3,5]),
}
def P3(V, w):
    Vm = [(mp.mpf(sp.N(x, 80)), mp.mpf(sp.N(y, 80))) for x, y in V]
    H = direct.setup(Vm, [mp.mpf(sp.N(x, 80)) for x in w])
    js = [direct.edge_jet_fast(H, q) for q in range(3)]
    return [2*sum(j[k] for j in js) for k in range(4)]
for name, (V, w) in cases.items():
    eds = formula.edges_data(V, w)
    P = P3(V, w)
    F = 2*sum(formula.E_formula(ed) for ed in eds)
    print('==', name, ' ratios r =', [ed['r'] for ed in eds])
    print("   P'''(1) direct =", mp.nstr(P[3], 45), '  table - direct =', mp.nstr(F - P[3], 3))
    if 'affine' in name:
        # hexagon law: P(1+d) = 2^{d/(1+d)} [P(1) + pi^2 E/24 d^2 - E/8 (z3 + 2pi^2/3) d^3] + O(d^4), E = sum_e h_e^2/L_e
        Ecal = 2*sum((ed['L2']*ed['V02'] - ed['q0']**2)/ed['L2']/sp.sqrt(ed['L2']) for ed in eds)
        d = sp.symbols('d')
        law = sp.exp(sp.log(2)*d/(1+d))*(sum(2*sp.sqrt(ed['L2']) for ed in eds) + sp.pi**2*Ecal/24*d**2 - Ecal/8*(sp.zeta(3)+2*sp.pi**2/3)*d**3)
        c3 = sp.series(law, d, 0, 4).removeO().coeff(d, 3)*6
        print("   hexagon law P'''(1) =", mp.nstr(mp.mpf(sp.N(c3, 70)), 45), '  law - direct =', mp.nstr(mp.mpf(sp.N(c3, 70)) - P[3], 3))
