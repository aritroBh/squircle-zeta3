import mpmath as mp, json
mp.mp.dps = 120
from fractions import Fraction
from math import factorial, comb

# exact H_n
N=180
KK=N+2
a=[Fraction(0)]*KK; a[0]=Fraction(1)
for k in range(1,KK): a[k]=Fraction((-2)**k,2*factorial(k))
s=[Fraction(0)]*KK; s[0]=Fraction(1)
for k in range(1,KK): s[k]=(a[k]-sum(s[i]*s[k-i] for i in range(1,k)))/2
Hmp=[mp.mpf((2*s[n]).numerator)/mp.mpf((2*s[n]).denominator) for n in range(N+1)]
f450=[mp.mpf(v) for v in json.load(open('../D11b2_stokes_jump/step2_num450_N180.json'))['f450']]

# r_n = f_n/n! - H_n
r = [f450[n]/mp.factorial(n) - Hmp[n] for n in range(N+1)]
# root test on r_n: radius should be pi (2*rho)
import math
print("root test |r_n|^{1/n} vs 1/pi =", mp.nstr(1/mp.pi,10))
for n in [60, 100, 140, 180]:
    print(f" n={n}: |r_n|^(1/n) = {mp.nstr(abs(r[n])**(mp.mpf(1)/n),10)}")

# beta_n = sum_{j>=0} r_{n-2j} (2 pi)^{n-2j} binom(n-j-1, j), n>=1 ; beta_0 = r_0
beta = [mp.mpf(0)]*(N+1)
beta[0] = r[0]
twopi = 2*mp.pi
for n in range(1, N+1):
    tot = mp.mpf(0)
    for j in range(0, (n-1)//2 + 1):
        k = n-2*j
        if k < 1: continue
        tot += r[k]*(twopi**k)*comb(n-j-1, j)
    beta[n] = tot
# sanity: evaluate R(zeta) on the disk vs direct series
def R_direct(z):
    return sum(r[n]*z**n for n in range(N+1))
def w_of_z(z):
    return z/(mp.pi + mp.sqrt(mp.pi**2 + z**2))
def R_map(z):
    w = w_of_z(z)
    return sum(beta[n]*w**n for n in range(N+1))
for z in [0.5, 1.0, 2.0, 3.0]:
    print(f"z={z}: direct={mp.nstr(R_direct(mp.mpf(z)),16)}  map={mp.nstr(R_map(mp.mpf(z)),16)}  diff={mp.nstr(abs(R_direct(mp.mpf(z))-R_map(mp.mpf(z))),4)}")
with open('progress.log','a') as f:
    f.write("step_16: H-split conformal map built (zeta=2 pi w/(1-w^2) exact); disk sanity vs direct series\n")
