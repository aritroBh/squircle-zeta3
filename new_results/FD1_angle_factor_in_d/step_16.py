lib = r'''
"""fd1_lib.py -- FD1 core: F_d(c) via the ALPHA2 leading law + corner_series (no kernels)."""
import sys, os
sys.path.insert(0, '<project>/research/08_theory/FD1_angle_factor_in_d')
os.chdir('<project>/research/08_theory/FD1_angle_factor_in_d')
from mpmath import mp, mpf, nstr, pi, sin, cos, acos, gamma, hyper, factorial
import corner_series

A_ZIG = [1,1,1,2,5,16,61,272,1385,7936,353792,22368256,1903757312,209865342976,
         29088885112832, 4951498053124096]

def eps(m): return (-1)**((m+1)//2)

def Eafun_factory(ca):
    def Ea(s):
        c2x = (1+ca)/2
        return gamma(1-s)*gamma(2+s)*hyper([1+s, -s],[2], c2x)/(1-ca)
    return Ea

def ew_coeffs(ca, nmax):
    """[y^w] E_a(y) for w=0..nmax via series of hypergeometric + gamma (exact formal)."""
    Ea = Eafun_factory(ca)
    from mpmath import diffs
    return [list(diffs(Ea, 0, n))[-1]/factorial(n) for n in range(nmax+1)]

def F_from_law(a, d, ew=None):
    """Numeric F_d(c) at angle a from the leading law: c(w;a) with w=d+1."""
    w = d+1
    ca = cos(a)
    if ew is None: ew = ew_coeffs(ca, w)
    M = corner_series.corner_coeffs(a, w-1)   # M_0..M_{w-1}; c(w)=M_{w-1}
    cw = M[w-1]
    trig = (1-ca)/2 if w % 2 == 0 else (1+ca)/2
    pref = eps(w-2)*A_ZIG[w-1]*trig*ew[w]/factorial(w-2)
    return cw/pref
'''
fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
open(fd1+'/fd1_lib.py','w').write(lib)

# quick validation at one angle, low precision first: a with c=cos(a/2)=1/2 -> a/2 = pi/3 -> a=2pi/3
from mpmath import mp, mpf, nstr, pi, acos
mp.dps=50
import sys; sys.path.insert(0,fd1)
from fd1_lib import F_from_law
import time
c = mpf(1)/2
a = 2*acos(c)   # c = cos(a/2)
for d in range(2,6):
    t0=time.time()
    F = F_from_law(a, d)
    print(d, nstr(F,30), f'{time.time()-t0:.1f}s')
