import sys, os
workdir = '<project>/research/08_theory/TOW2R_review'
sys.path.insert(0, os.path.join(workdir, 'tow2_copy'))
os.chdir(os.path.join(workdir, 'tow2_copy'))
from fractions import Fraction as Fr
import genintegrand as GI
import sympy as sp
import tow2lib as T

def qi_mul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def qi_div(x,y):
    n = y[0]**2+y[1]**2
    return ((x[0]*y[0]+x[1]*y[1])/n, (x[1]*y[0]-x[0]*y[1])/n)
TSTAR = (Fr(1,2), Fr(1,2))
def peval(P, x):
    out = (Fr(0),Fr(0)); xp = (Fr(1),Fr(0))
    for c in P:
        out = (out[0]+c*xp[0], out[1]+c*xp[1])
        xp = qi_mul(xp, x)
    return out
def pderiv(P):
    return [c*(k+1) for k,c in enumerate(P[1:])] if len(P)>1 else [Fr(0)]

d_sym, m_sym = T.d, T.m
bad = []; tested = 0
for dd in range(4, 13):
    for mm in range((dd+1)//2, dd+1):
        S = GI.S_poly(dd, mm)
        cur = list(S); ders = [list(cur)]
        for l in range(1, 5):
            cur = pderiv(cur); ders.append(list(cur))
        Sv = peval(ders[0], TSTAR)
        for l in (3,4):
            val = qi_div(peval(ders[l], TSTAR), Sv)
            cl = sp.simplify(T.R[l]).subs({d_sym: dd, m_sym: mm})
            re, im = sp.re(cl), sp.im(cl)
            ok = (val[0] == Fr(re.numerator, re.denominator)) and (val[1] == Fr(im.numerator, im.denominator))
            tested += 1
            if not ok: bad.append((dd,mm,l))
print(f'S^l/S (l=3,4) direct-from-definition vs N_l/Pi_l: {tested} cells, {len(bad)} bad')
if bad: print(bad[:8])
with open(os.path.join(workdir,'progress.log'),'a') as f:
    f.write(f'step_10: independent S3/S4 (R_3,R_4) check from definition: {tested} cells d=4..12, {len(bad)} bad\n')
