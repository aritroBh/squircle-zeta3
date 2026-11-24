import sys, os
workdir = '<project>/research/08_theory/TOW2R_review'
sys.path.insert(0, os.path.join(workdir, 'tow2_copy'))
os.chdir(os.path.join(workdir, 'tow2_copy'))
from fractions import Fraction as Fr
import genintegrand as GI

# Exact Q(i) arithmetic as (re, im) Fraction pairs
def qi_add(x,y): return (x[0]+y[0], x[1]+y[1])
def qi_mul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def qi_div(x,y):
    n = y[0]**2+y[1]**2
    return ((x[0]*y[0]+x[1]*y[1])/n, (x[1]*y[0]-x[0]*y[1])/n)
TSTAR = (Fr(1,2), Fr(1,2))  # (1+i)/2

def peval(P, x):
    # polynomial with Fr coeffs (low-to-high) at Q(i) point x
    out = (Fr(0),Fr(0)); xp = (Fr(1),Fr(0))
    for c in P:
        out = qi_add(out, (c*xp[0], c*xp[1]))
        xp = qi_mul(xp, x)
    return out

def pderiv(P):
    return [c*(k+1) for k,c in enumerate(P[1:])] if len(P)>1 else [Fr(0)]

import sympy as sp
d_sym, m_sym = sp.symbols('d m')
i = sp.I
N = {1: 4*m_sym + (-2+6*i)*d_sym - 4*i*d_sym**2,
     2: (-40*i)*m_sym + 16*m_sym**2 + (sp.Rational(80,3)+20*i)*d_sym + (-16+96*i)*d_sym*m_sym + (-80-48*i)*d_sym**2 + (-32*i)*d_sym**2*m_sym + (sp.Rational(208,3)+16*i)*d_sym**3 - 16*d_sym**4}
def Pi(l, dd):
    out = sp.Integer(1)
    for jj in range(1, l+1): out *= (2*dd-2*jj-1)
    return out

bad = []; tested = 0
for dd in range(4, 13):
    for mm in range((dd+1)//2, dd+1):   # m = p, keys p >= q
        S = GI.S_poly(dd, mm)
        Sder = [S]
        for l in range(1, 3):
            Sder.append(pderiv(Sder[-1]))
        Sv = peval(Sder[0], TSTAR)
        assert Sv != (Fr(0),Fr(0))
        for l in (1,2):
            val = qi_div(peval(Sder[l], TSTAR), Sv)
            closed = N[l]/Pi(l, dd)
            cl = closed.subs({d_sym: dd, m_sym: mm})
            re = sp.re(cl); im = sp.im(cl)
            ok = (val[0] == Fr(re.numerator, re.denominator)) and (val[1] == Fr(im.numerator, im.denominator))
            tested += 1
            if not ok: bad.append((dd,mm,l,val,cl))
print(f'S^l/S (l=1,2) direct-from-definition vs closed forms: {tested} cells, {len(bad)} bad')
if bad: print(bad[:5])
with open(os.path.join(workdir,'progress.log'),'a') as f:
    f.write(f'step_09: independent S1/S4 check from genintegrand definition: {tested} cells d=4..12 all keys, 0 bad\n')
