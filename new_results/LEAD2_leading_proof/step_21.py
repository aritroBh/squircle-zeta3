import os, time
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
import mpmath
mpmath.mp.dps = 60
from mpmath import mpf, sqrt, log, pi, zeta, gamma, sin, quad
import sys; sys.path.insert(0,'.')
from lead2_lib import Efun
from math import comb

print("=== R1: Beta substitution  E(z) = int_0^1 (w^{1+z}+w^{-z})/(1+w^2)^{3/2} dw ===")
for z in [mpf('0.3'), mpf('-0.4'), mpf('0.73'), mpf('0.2j')]:
    lhs = Efun(z)
    rhs = quad(lambda w: (w**(1+z)+w**(-z))/(1+w**2)**mpf('1.5'), [0,1])
    print(f"  z={mpmath.nstr(z,8)}: diff={mpmath.nstr(lhs-rhs,4)}")

print("=== R2: IBP  E(z)-E(-z) = z (F(z)+F(-z)),  F(z)=int_0^1 w^{-z} h dw ===")
def hfun(w): return ((1+w)/sqrt(1+w**2)-1)/w
def F(z): return quad(lambda w: w**(-z)*hfun(w), [0,1])
for z in [mpf('0.3'), mpf('-0.25'), mpf('0.6'), mpf('0.31j')]:
    lhs = Efun(z)-Efun(-z); rhs = z*(F(z)+F(-z))
    print(f"  z={mpmath.nstr(z,8)}: diff={mpmath.nstr(lhs-rhs,4)}")

print("=== R3: F(y) partial fractions: check n>=0 vs n>=1 in second sum ===")
def cn(n): return mpf(comb(2*n,n))/4**n
def F_series_mine(y):
    return sum((-1)**n*cn(n)/(2*n-y) for n in range(1,400)) + \
           sum((-1)**n*cn(n)/(2*n+1-y) for n in range(0,400))
def F_series_lead1(y):
    return sum((-1)**n*cn(n)*(1/(2*n-y)+1/(2*n+1-y)) for n in range(1,400))
y = mpf('0.3')
fq = F(y)
print("  F(0.3) quad   =", mpmath.nstr(fq,25))
print("  mine   (n>=0) :", mpmath.nstr(F_series_mine(y)-fq,6))
print("  lead1  (n>=1) :", mpmath.nstr(F_series_lead1(y)-fq,6))

print("=== R4: involution recursion a~_{p,q} = a~_{q+1,p-1} via A-integrals ===")
from mpmath import exp as mexp
def A(p,q):
    return (-1)**p*quad(lambda r: r**p*log(1-mexp(-r))**q/(mexp(r)-1), [0,1,mpmath.inf])
import math
def at(p,q): return A(p,q)/(mpmath.factorial(p)*mpmath.factorial(q))
bad = 0
for p in range(1,6):
    for q in range(0,5):
        r = at(p,q) - at(q+1,p-1)
        if abs(r) > mpf(10)**(-40): bad += 1; print("  BAD", p,q,r)
print("  recursion checked p=1..5,q=0..4; violations:", bad)
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: referee R1 Beta-sub, R2 IBP, R3 F-series n>=0 (LEAD1 formula has n>=1 erratum), R4 recursion\n")
