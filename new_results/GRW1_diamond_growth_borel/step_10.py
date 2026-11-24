
import mpmath as mp
mp.mp.dps = 60
from fractions import Fraction
from math import factorial
import numpy as np

N = 240
def Hseries_rhombus(q, N):
    a = [(Fraction(1,factorial(n)) if n%2==0 else q/Fraction(factorial(n))) for n in range(N+1)]
    g = [Fraction(0)]*(N+1); g[0]=Fraction(1)
    for k in range(1, N+1):
        g[k] = (a[k] - sum(g[i]*g[k-i] for i in range(1,k)))/2
    H = [Fraction(0)]*(N+1)
    for n in range(N+1):
        H[n] = sum((g[k] + ((-1)**k)*g[k]) * Fraction((-1)**(n-k), 2**(n-k)*factorial(n-k)) for k in range(0, n+1))
    return H

def darb_terms(zlist, clist, n, sgn):
    # -(1/(2 sqrt pi)) n^-3/2 sum_j c_j (-z_j)^{1/2} z_j^{-n} (1 - 3/(8n) + sgn*3 z_j/(4n))
    tot = mp.mpc(0)
    for z, c in zip(zlist, clist):
        tot += c*mp.sqrt(-z)*z**(-n)*(1 - mp.mpf(3)/(8*n) + sgn*3*z/(4*n))
    return (-1/(2*mp.sqrt(mp.pi)))*n**(-mp.mpf('1.5'))*tot

# ---- square ----
H0 = Hseries_rhombus(Fraction(0), N)
H0n = [mp.mpf(v.numerator)/mp.mpf(v.denominator) for v in H0]
z0, z1 = mp.mpc(0, mp.pi/2), mp.mpc(0, -mp.pi/2)
for sgn in [+1, -1]:
    errs = []
    for n in range(120, 241, 4):
        law2 = darb_terms([z0, z1], [2, 2], n, sgn).real
        errs.append(abs(1 - H0n[n]/law2))
    print(f"square sgn={sgn:+d}: |1-H_n/law2| at n=120: {mp.nstr(errs[0],4)}, n=240: {mp.nstr(errs[-1],4)}, ratio~n^2: {mp.nstr(errs[-1]*240**2,6)}")

# ---- rhombus q=-3/5 ----
Hr = Hseries_rhombus(Fraction(-3,5), N)
Hrn = [mp.mpf(v.numerator)/mp.mpf(v.denominator) for v in Hr]
A, B = 1, 2; r = mp.sqrt(5); L = mp.log(2); rho = mp.sqrt(L**2+mp.pi**2/4)
th0 = mp.atan2(mp.pi/2, L)
zlist = [L + 1j*mp.pi/2, -L + 1j*mp.pi/2, L - 1j*mp.pi/2, -L - 1j*mp.pi/2]
clist = [A*mp.sqrt(2)/r, B*mp.sqrt(2)/r, A*mp.sqrt(2)/r, B*mp.sqrt(2)/r]
def bracket(n):
    return B*((-1)**n)*mp.cos((n-mp.mpf('0.5'))*th0) - A*mp.sin((n-mp.mpf('0.5'))*th0)
for sgn in [+1, -1]:
    errs = []
    for n in range(120, 241):
        if abs(bracket(n)) < 1.0: continue
        law2 = darb_terms(zlist, clist, n, sgn).real
        errs.append((n, abs(1 - Hrn[n]/law2)))
    print(f"rhombus sgn={sgn:+d}: worst |1-ratio| over safe n: {mp.nstr(max(e for _,e in errs),5)}; at n={errs[-1][0]}: {mp.nstr(errs[-1][1],4)} x n^2={mp.nstr(errs[-1][1]*errs[-1][0]**2,5)}")
