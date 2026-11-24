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

def ratioR(d, m, l):
    return S_val(d,m,l)/S_val(d,m,0)

dd, mm = sp.symbols('d m')
def fit(l, dmax=13):
    mons = [(a,b) for a in range(0,2*l+1) for b in range(0,2*l+1-a)]
    Arows = []; br_re = []; br_im = []; pts = []
    for d in range(l+1, dmax+1):
        for m in range(0, d+1):
            R = ratioR(d,m,l)
            prod = 1
            for j in range(1,l+1): prod *= (2*d-2*j-1)
            br_re.append(sp.Rational(R.p0)*prod); br_im.append(sp.Rational(R.p2)*prod)
            Arows.append([sp.Integer(d)**a * sp.Integer(m)**b for (a,b) in mons]); pts.append((d,m))
    Amat = sp.Matrix(Arows); n = len(mons); idxs = list(range(len(pts)))
    for attempt in range(200):
        sub = idxs[attempt*n:(attempt+1)*n]
        if len(sub) < n: break
        M = Amat[sub, :]
        if M.det() != 0:
            sol_re = M.solve(sp.Matrix([br_re[t] for t in sub]))
            sol_im = M.solve(sp.Matrix([br_im[t] for t in sub]))
            ok = all(Amat[t,:]*sol_re == sp.Matrix([br_re[t]]) and Amat[t,:]*sol_im == sp.Matrix([br_im[t]]) for t in idxs)
            if ok:
                poly_re = sp.expand(sum(sol_re[q]*dd**mons[q][0]*mm**mons[q][1] for q in range(n)))
                poly_im = sp.expand(sum(sol_im[q]*dd**mons[q][0]*mm**mons[q][1] for q in range(n)))
                return poly_re, poly_im, len(pts)
    return None, None, 0

fits = {}
for l in range(1,4):
    pr, pi, npts = fit(l)
    fits[l] = (pr, pi)
    print(f'=== l={l}: numerator fit over {npts} points (poly in (d,m) of degree {2*l})')
    print('  RE:', sp.factor(pr))
    print('  IM:', sp.factor(pi))
import pickle
pickle.dump({l:(str(fits[l][0]),str(fits[l][1])) for l in fits}, open('Sderiv_fits.pkl','wb'))
