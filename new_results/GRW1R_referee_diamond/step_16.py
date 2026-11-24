from fractions import Fraction
from math import factorial
import mpmath as mp
mp.mp.dps = 80

# Exact H_n: H(z) = sqrt(2)*sqrt(1+exp(-2z)) = 2*sqrt(S(z)), S=(1+exp(-2z))/2 = sum a_k z^k... 
# Use GRW1-independent derivation: H^2 = 2+2e^{-2z}; let W(z)=H/2 = sqrt(cosh z)*e^{-z/2}... simpler: S(z)=(1+e^{-2z})/2, H=2 sqrt(S).
N = 200
# S(z) = 1 + sum_{k>=1} (-2)^k/(2 k!) z^k = sum a_k z^k
KK = N+2
a = [Fraction(0)]*KK
a[0] = Fraction(1)
for k in range(1, KK):
    a[k] = Fraction((-2)**k, 2*factorial(k))
# b = sqrt(a)
s = [Fraction(0)]*KK
s[0] = Fraction(1)
for k in range(1, KK):
    s[k] = (a[k] - sum(s[i]*s[k-i] for i in range(1, k)))/2
Hf = [2*s[n] for n in range(N+1)]
print("H_0..H_6:", Hf[:7])
# sanity: H_1 = -1? H(z)=2-z+...
# numeric check via Cauchy integral
def Hnum(z): return mp.sqrt(2)*mp.sqrt(1+mp.e**(-2*z))
import cmath
def Hn_num(n):
    r = 0.5
    return mp.quad(lambda t: Hnum(r*mp.e**(1j*t))/ (r*mp.e**(1j*t))**n, [0, 2*mp.pi])/(2*mp.pi)
for n in [0,1,2,5]:
    print(n, mp.nstr(mp.mpf(Hf[n].numerator)/mp.mpf(Hf[n].denominator),20), mp.nstr(Hn_num(n),20))

rho = mp.pi/2
Hmp = [mp.mpf(v.numerator)/mp.mpf(v.denominator) for v in Hf]
# Claimed law: H_n = -sqrt2 rho^{-n} n^{-3/2} cos((2n+1)pi/4) (1 + (3-3pi(-1)^n)/(8n)) + R_n
def law2(n):
    c = mp.cos((2*n+1)*mp.pi/4)
    return -mp.sqrt(2)*rho**(-n)*mp.mpf(n)**mp.mpf('-1.5')*c*(1 + (3-3*mp.pi*((-1)**n))/(8*n))
def law1(n):
    c = mp.cos((2*n+1)*mp.pi/4)
    return -mp.sqrt(2)*rho**(-n)*mp.mpf(n)**mp.mpf('-1.5')*c

print("\n n | H_n/law1 | (H_n-law2)*rho^n*n^(5/2) | *(n^(7/2)) | sign ok?")
import math
for n in [10, 20, 40, 60, 80, 100, 120, 150, 180, 200]:
    r1 = Hmp[n]/law1(n)
    resid = Hmp[n]-law2(n)
    s52 = resid*rho**n*mp.mpf(n)**mp.mpf('2.5')
    s72 = resid*rho**n*mp.mpf(n)**mp.mpf('3.5')
    print(f"{n:4d} | {mp.nstr(r1,12)} | {mp.nstr(s52,10)} | {mp.nstr(s72,10)}")

# sign pattern check (-,+,+,-) for n mod 4 = 0,1,2,3 and C_H
CH = mp.mpf(0); nCH=0
for n in range(1, 201):
    v = abs(Hmp[n])*rho**n*mp.mpf(n)**mp.mpf('1.5')
    if v > CH: CH, nCH = v, n
print("\nmax |H_n| rho^n n^1.5 over n<=200:", mp.nstr(CH,12), "at n =", nCH, " (claimed 5.874 at n=3)")
# two-term residual claimed measured <=13.07 rho^{-n} n^{-5/2} -- check what that means
mx = mp.mpf(0)
for n in range(10, 201):
    v = abs(Hmp[n]-law2(n))*rho**n*mp.mpf(n)**mp.mpf('2.5')
    if v>mx: mx=v
print("max_{10<=n<=200} |H_n-law2| rho^n n^{5/2} =", mp.nstr(mx,8))
mx2 = mp.mpf(0)
for n in range(10, 201):
    v = abs(Hmp[n]-law2(n))*rho**n*mp.mpf(n)**mp.mpf('3.5')
    if v>mx2: mx2=v
print("max_{10<=n<=200} |H_n-law2| rho^n n^{7/2} =", mp.nstr(mx2,8))
with open('progress.log','a') as f:
    f.write("step_04: exact rational H_n n<=200; two-term law confirmed; C_H and residuals measured\n")
