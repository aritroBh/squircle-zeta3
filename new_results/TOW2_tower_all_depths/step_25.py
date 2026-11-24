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

l = 1
mons = [(a,b) for a in range(0,2*l+1) for b in range(0,2*l+1-a)]
print('mons:', mons)
rows = []; vals = []
for d in range(2, 8):
    for m in range(0, d+1):
        R = S_val(d,m,l)/S_val(d,m,0)
        prod = (2*d-3)
        rows.append([sp.Integer(d)**a*sp.Integer(m)**b for (a,b) in mons])
        vals.append(sp.Rational(R.p0)*prod + sp.I*sp.Rational(R.p2)*prod)
Amat = sp.Matrix(rows)
print('rank of Amat:', Amat.rank(), 'rows:', Amat.shape)
# solve complex directly
n = len(mons)
M = Amat[:n,:]
print('det:', M.det())
sol = M.solve(sp.Matrix(vals[:n]))
resid = Amat*sol - sp.Matrix(vals)
print('max residual:', max(abs(complex(r)) for r in resid))
dd, mm = sp.symbols('d m')
poly = sp.expand(sum(sol[q]*dd**mons[q][0]*mm**mons[q][1] for q in range(n)))
print('N_1(d,m) =', sp.factor(poly))
print('compare S3: ', sp.expand((2*(dd-2*mm)-2*dd*sp.I)*(2*dd-3) - 4*(dd-1)*(dd-2*mm)))
