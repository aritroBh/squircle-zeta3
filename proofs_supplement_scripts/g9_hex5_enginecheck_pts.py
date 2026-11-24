# g9_hex5_enginecheck_pts.py -- machine step M0 at rational points: the d^k integrand, k <= 5, of a hexagon edge by the engine of Theorem 25
# (g4_edge.py) against the table E5_polys.pkl of the notes, all coefficients of t^n X^a Y^b W^c, at three rational points of the edge data
# (all arithmetic exact).  The symbolic version is g9_hex5_enginecheck.py.
import os, sys, pickle, time, sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g4_edge import edge_E
A0, A1, B, kap = sp.symbols('A0 A1 B kap'); L2, q0, V02 = sp.symbols('L2 q0 V02'); lA0, lA1, lB = sp.symbols('lA0 lA1 lB')
t, X, Y, W = sp.symbols('t X Y W')
R = sp.Rational
E5FILE = os.environ.get('G9E5', os.path.join(os.environ.get('G9WORK', '.'), 'E5_polys.pkl'))
notes = pickle.load(open(E5FILE, 'rb'))
pts = [dict(B=R(1, 3), kap=R(2, 5), q0=R(-2, 7), V02=R(3, 5), Lq=R(3, 2), lA0=R(3, 11), lA1=R(-5, 13), lB=R(7, 9)),
       dict(B=R(2, 7), kap=R(-3, 4), q0=R(5, 3), V02=R(7, 2), Lq=R(5, 4), lA0=R(-1, 3), lA1=R(2, 5), lB=R(-4, 7)),
       dict(B=R(3, 5), kap=R(1, 9), q0=R(-1, 2), V02=R(9, 4), Lq=R(2, 3), lA0=R(5, 6), lA1=R(1, 7), lB=R(3, 2))]
t0 = time.time()
for P in pts:
    Bv, kv = P['B'], P['kap']; Lv = P['Lq']
    a0, a1, a2 = (1-Bv-Bv*kv)*t, (1-Bv)*(1-t), Bv*(1+kv*t)
    terms = [(a0, 1-Bv-Bv*kv, X+P['lA0']), (a1, -(1-Bv), Y+P['lA1']), (a2, Bv*kv, W+P['lB'])]
    E = edge_E(terms, Lv, P['q0'] + Lv**2*t, P['V02'] + 2*P['q0']*t + Lv**2*t**2, 5)
    sub = {A1: 1-Bv, A0: 1-Bv-Bv*kv, B: Bv, kap: kv, q0: P['q0'], V02: P['V02'], L2: Lv**2, lA0: P['lA0'], lA1: P['lA1'], lB: P['lB']}
    bad = 0; tot = 0
    for k in range(6):
        mine = {m: c for m, c in sp.Poly(sp.expand(E[k]), t, X, Y, W).terms()}
        theirs = {}
        for m, c in notes[k].items():
            v = sp.sympify(c).subs(sub)
            if v != 0: theirs[m] = v
        for m in set(mine) | set(theirs):
            tot += 1
            if sp.simplify(mine.get(m, 0) - theirs.get(m, 0)) != 0: bad += 1
    print('point B=%s kappa=%s: %d coefficients (orders 0..5) compared, differing: %d  (%.0f s)' % (Bv, kv, tot, bad, time.time()-t0))
