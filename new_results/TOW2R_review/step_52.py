import sys, os
workdir = '<project>/research/08_theory/TOW2R_review'
sys.path.insert(0, os.path.join(workdir, 'tow2_copy'))
os.chdir(os.path.join(workdir, 'tow2_copy'))
import mpmath as mp
mp.mp.dps = 50
from math import comb, factorial
import genintegrand as GI

def t_of_v(v): return (v**2+2*v-1)/(4*v)
def q2_of_t(t): return t**2+(1-t)**2
def S2(n,k):
    if n<0 or k<0: return mp.mpf(0)
    if n==0 and k==0: return mp.mpf(1)
    if k==0 or n==0 or k>n: return mp.mpf(0)
    return k*S2(n-1,k)+S2(n-1,k-1)
def M_k_sqrt_at_v(v, d, m, sign):
    t = t_of_v(v); q2 = q2_of_t(t); q = mp.sqrt(q2)
    alpha=(1-t)**2/q2; beta=t**2/q2
    tot = mp.mpf(0)
    for j in range(1, d+1):
        bj = mp.binomial(mp.mpf('0.5'), j)
        for s in range(0, j+1):
            if m < s or d-m < j-s: continue
            cx = factorial(s)*S2(m,s)*2**m/factorial(m)
            cy = factorial(j-s)*S2(d-m,j-s)*2**(d-m)/factorial(d-m)
            tot += bj*comb(j,s)*alpha**s*beta**(j-s)*cx*cy
    Qp = q/(t*(1-t))
    dtdv = (v**2+1)/(4*v**2)
    return sign*mp.mpf('0.5')*Qp*tot*dtdv
def Sval(S, t):
    return sum(mp.mpf(c.numerator)/c.denominator * t**k for k,c in enumerate(S))
def formula_coeff(v, d, m):
    t = t_of_v(v); q2 = q2_of_t(t); q = mp.sqrt(q2)
    S = GI.S_poly(d, m)
    dtdv = (v**2+1)/(4*v**2)
    return -q*Sval(S,t)/(2*t*(1-t)*q2**d)*dtdv

v = mp.mpf('0.7')+mp.mpf('0.13')*1j
for (d,m) in [(4,2),(4,4),(5,3),(6,1)]:
    a = M_k_sqrt_at_v(v,d,m,+1); b = formula_coeff(v,d,m)
    print((d,m), 'ratio a/b =', mp.nstr(a/b, 6))
