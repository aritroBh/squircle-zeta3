# g9_hex5_enginecheck.py -- machine step "E5": the d^5 integrand of a hexagon edge, by the engine of Theorem 25 (g4_edge.py, K = 5),
# compared monomial by monomial with the table E5_polys.pkl rebuilt from the notes' engine (g9_hex5_s1_engine.py).
# Usage: G9WORK=<dir with E5_polys.pkl> python g9_hex5_enginecheck.py   (or G9E5=<path of E5_polys.pkl>)
import os, sys, pickle, time, sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g4_edge import edge_E
WD = os.environ['G9WORK']
t0 = time.time()
A0, A1, B, kap = sp.symbols('A0 A1 B kap'); L2, q0, V02 = sp.symbols('L2 q0 V02'); lA0, lA1, lB = sp.symbols('lA0 lA1 lB')
t, X, Y, W = sp.symbols('t X Y W')
Lp = sp.Symbol('L', positive=True)     # L; the notes write sqrt(L2)
Lsym = Lp
rel = {A1: 1 - B, A0: 1 - B - B*kap}
a0, a1, a2 = (1-B-B*kap)*t, (1-B)*(1-t), B*(1+kap*t)
terms = [(a0, 1-B-B*kap, X+lA0), (a1, -(1-B), Y+lA1), (a2, B*kap, W+lB)]
E = edge_E(terms, Lsym, q0 + Lp**2*t, V02 + 2*q0*t + Lp**2*t**2, 5)
mine = {}
for k in range(6):
    P = sp.Poly(sp.expand(E[k]), t, X, Y, W)
    mine[k] = {m: c for m, c in P.terms()}
print('g4 engine, K=5: monomials per order', {k: len(v) for k, v in mine.items()}, '(%.0f s)' % (time.time()-t0))
E5FILE = os.environ.get('G9E5', f'{WD}/E5_polys.pkl')
notes = pickle.load(open(E5FILE, 'rb'))
for k in range(6):
    keys = set(mine[k]) | set(notes[k]); bad = 0
    for m in keys:
        a = mine[k].get(m, 0); b = sp.expand(sp.sympify(notes[k].get(m, 0)).subs(rel).subs(L2, Lp**2))
        if sp.expand(a - b) != 0: bad += 1
    print('order', k, ': monomials t^n X^a Y^b W^c mine', len(mine[k]), ' notes', len(notes[k]), ' differing coefficients', bad)
print('done (%.0f s)' % (time.time()-t0))
