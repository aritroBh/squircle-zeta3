import os, time, math
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
import mpmath
mpmath.mp.dps = 60
from mpmath import mpf, mpc, sqrt, log, pi, zeta, gamma, sin, quad, exp as mexp
import sys; sys.path.insert(0,'.')
from lead2_lib import Efun
from math import comb

print("=== R1: Beta substitution  E(z) = int_0^1 (w^{1+z}+w^{-z})/(1+w^2)^{3/2} dw ===")
for z in [mpf('0.3'), mpf('-0.4'), mpf('0.73'), mpc(0, '0.2')]:
    lhs = Efun(z)
    rhs = quad(lambda w: (w**(1+z)+w**(-z))/(1+w**2)**mpf('1.5'), [0,1])
    print(f"  z={mpmath.nstr(z,8)}: diff={mpmath.nstr(abs(lhs-rhs),4)}")

print("=== R2: IBP  E(z)-E(-z) = z (F(z)+F(-z)) ===")
def hfun(w): return ((1+w)/sqrt(1+w**2)-1)/w
def F(z): return quad(lambda w: w**(-z)*hfun(w), [0,1])
for z in [mpf('0.3'), mpf('-0.25'), mpf('0.6'), mpc(0,'0.31')]:
    lhs = Efun(z)-Efun(-z); rhs = z*(F(z)+F(-z))
    print(f"  z={mpmath.nstr(z,8)}: diff={mpmath.nstr(abs(lhs-rhs),4)}")

print("=== R3: F(y) partial fractions: n>=0 (mine) vs n>=1 (LEAD1 genfun.md print) ===")
def cn(n): return mpf(comb(2*n,n))/4**n
def F_series_mine(y):
    return sum((-1)**n*cn(n)/(2*n-y) for n in range(1,600)) + \
           sum((-1)**n*cn(n)/(2*n+1-y) for n in range(0,600))
def F_series_lead1(y):
    return sum((-1)**n*cn(n)*(1/(2*n-y)+1/(2*n+1-y)) for n in range(1,600))
y = mpf('0.3'); fq = F(y)
print("  mine   (n>=0) residual:", mpmath.nstr(F_series_mine(y)-fq,6))
print("  lead1  (n>=1) residual:", mpmath.nstr(F_series_lead1(y)-fq,6))

print("=== R4: involution recursion a~_{p,q} = a~_{q+1,p-1} ===")
def A(p,q):
    return (-1)**p*quad(lambda r: r**p*log(1-mexp(-r))**q/(mexp(r)-1), [0,1,mpmath.inf])
def at(p,q): return A(p,q)/(mpmath.factorial(p)*mpmath.factorial(q))
worst = mpf(0)
for p in range(1,6):
    for q in range(0,5):
        worst = max(worst, abs(at(p,q)-at(q+1,p-1)))
print("  max |residual| p=1..5,q=0..4:", mpmath.nstr(worst,4))
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: referee R1,R2,R3,R4 (Beta-sub, IBP, F-series erratum n>=0, recursion)\n")
