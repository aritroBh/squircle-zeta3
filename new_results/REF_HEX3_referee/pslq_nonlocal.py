"""Independent test of prop:nonlocal: per-edge P'''(1) contribution from the referee's direct jet integration,
then an integer relation on [E*L, z3, pi^2 log2, log^3 2, pi^2, log^2 2, log2, 1]; the zeta(3) coefficient
minus Z_e (printed prop:edgelaw) is compared with the printed N_e = Phi_A z3[A(1-r)] + Phi_L z3[Li3(1-r)]."""
import sympy as sp, mpmath as mp, direct, formula
from run_tests import frame
mp.mp.dps = 120
CASES = {
 'H2 paper example edge0 (r=2)': ([(1,0),(2,2),(0,1)], [1,3,1], 0),
 'H2 paper example edge1 (r=1/2)': ([(1,0),(2,2),(0,1)], [1,3,1], 1),
 'H2 paper example edge2 (r=1)': ([(1,0),(2,2),(0,1)], [1,3,1], 2),
 'H7 fresh edge0 (r=2)': (frame(2, sp.Rational(3,2), [[2,1],[1,3]]), [1,5,3], 0),
 'H8 fresh edge1 (r=1/2)': (frame(sp.Rational(3,2), 2, [[1,1],[-1,2]]), [3,5,1], 1),
}
l2 = mp.log(2)
for name, (V, w, q) in CASES.items():
    eds = formula.edges_data(V, w); ed = eds[q]
    print('==', name, 'r=', ed['r'], 'A0,A1,B=', ed['A0'], ed['A1'], ed['B'], 'A/w:', ed['A0']/ed['wA0'], ed['A1']/ed['wA1'], ed['B']/ed['wB'])
    Vm = [(mp.mpf(sp.N(x, 140)), mp.mpf(sp.N(y, 140))) for x, y in V]
    H = direct.setup(Vm, [mp.mpf(x) for x in w])
    E = direct.edge_jet_fast(H, q)[3]
    L = mp.sqrt(mp.mpf(sp.N(ed['L2'], 140)))
    X = E * L
    basis = [X, mp.zeta(3), mp.pi**2*l2, l2**3, mp.pi**2, l2**2, l2, 1]
    rel = mp.pslq(basis, maxcoeff=10**12, maxsteps=10**6)
    print('   relation', rel, ' residual', mp.nstr(mp.fsum(a*b for a, b in zip(rel, basis)), 5))
    cz3 = sp.Rational(-rel[1], rel[0]) / sp.sqrt(ed['L2'])
    Ze = formula.Ze_printed(ed); Ne = formula.Ne_printed(ed)
    print('   z3 coefficient (PSLQ):', sp.nsimplify(sp.radsimp(cz3)), ' Z_e (printed):', Ze, ' N_e (printed):', Ne,
          '  PSLQ - (Z_e+N_e) =', sp.simplify(cz3 - Ze - Ne))
