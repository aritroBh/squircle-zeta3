import sys, os
sys.dont_write_bytecode = True
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
from towlib import *
from genintegrand import S_poly
from fractions import Fraction as Fr
import sympy as sp

tstar = QSI(Fr(1,2),0,Fr(1,2),0)
def S_val(d, m, l):
    P = S_poly(d, m)
    out = ZERO
    for n, c in enumerate(P):
        if n < l: continue
        ff = 1
        for a in range(l): ff *= (n-a)
        out = out + QSI(Fr(c*ff)) * (tstar**(n-l))
    return out

dd, mm = sp.symbols('d m')
def fitR(l, dmax=13, deg=None):
    if deg is None: deg = 2*l
    mons = [(a,b) for a in range(0,deg+1) for b in range(0,deg+1-a)]
    rows = []; vals = []
    for d in range(l+1, dmax+1):
        for m in range(0, d+1):
            R = S_val(d,m,l)/S_val(d,m,0)
            prod = 1
            for j in range(1,l+1): prod *= (2*d-2*j-1)
            rows.append([sp.Integer(d)**a*sp.Integer(m)**b for (a,b) in mons])
            vals.append(sp.Rational(R.p0)*prod + sp.I*sp.Rational(R.p2)*prod)
    Amat = sp.Matrix(rows); b = sp.Matrix(vals)
    G = Amat.T*Amat
    if G.det() == 0:
        return None, 'singular normal eq'
    sol = G.solve(Amat.T*b)
    resid = Amat*sol - b
    maxres = max(abs(complex(r)) for r in resid)
    poly = sp.expand(sum(sol[q]*dd**mons[q][0]*mm**mons[q][1] for q in range(len(mons))))
    return poly, maxres

for l in range(1,5):
    poly, res = fitR(l)
    print(f'=== l={l}, residual {res}')
    print('   ', sp.factor(poly))
