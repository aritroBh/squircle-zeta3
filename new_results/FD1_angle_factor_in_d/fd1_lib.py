
"""fd1_lib.py -- FD1 core: F_d(c) via the ALPHA2 leading law + corner_series (no kernels)."""
import sys, os
sys.path.insert(0, '<project>/research/08_theory/FD1_angle_factor_in_d')
os.chdir('<project>/research/08_theory/FD1_angle_factor_in_d')
from mpmath import mp, mpf, nstr, pi, sin, cos, acos, gamma, hyper, factorial
import corner_series

A_ZIG = [1,1,1,2,5,16,61,272,1385,7936,50521,353792,2702765,199360981,
         19391512145, 2404879675441]

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

# ---- squircle top table (LEAD1 result.md sec.4 closed form), beta = D/(d-1)! ----
from fractions import Fraction as Fr
from math import factorial, comb
A_Z = [1,1,1,2,5,16,61,272,1385,7936,50521,353792,2702765,199360981,19391512145,2404879675441]
def eps_m(m): return (-1)**((m+1)//2)
def delta_d(d):  # [x^d](sech+tanh) = eps_{d-1} A_d/d!
    return Fr(eps_m(d-1)*A_Z[d], factorial(d))
def squircle_D(d, i, j, pole):
    """LEAD1 closed-form top table. pole in {'-a','a','-b','0'}."""
    dd = delta_d(d)
    lastkey = (i == (d+1)//2 and j == d//2)   # the 'last key'
    diag = (i == j)
    if pole == '-a':
        return Fr((-1)**(j+1))*comb(d,j)*dd
    if pole == 'a':
        if d % 2 == 0: return Fr(0)
        v = Fr((-1)**j)*comb(d,j)*dd
        return v/2 if lastkey else v
    if pole == '-b':
        if d % 2 == 0:
            if j < d/2: return Fr((-1)**j)*comb(d,j)*dd
            return Fr(0)
        # odd d: last key also carries -b with the NEGATIVE of the pole-0 value
        if lastkey: return -Fr((-1)**j)*comb(d,j)*dd/2
        return Fr(0)
    if pole == '0':
        if diag or (d%2==1 and lastkey):
            return Fr((-1)**j)*comb(d,j)*dd/2
        return Fr(0)
def squircle_beta(d, i, j, pole):
    return squircle_D(d,i,j,pole)/factorial(d-1)
