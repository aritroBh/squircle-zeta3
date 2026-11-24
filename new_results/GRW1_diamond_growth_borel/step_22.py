
import mpmath as mp
mp.mp.dps = 80
import json
from fractions import Fraction
from math import factorial as ffact

# ---- exact H_n ----
N = 300
KK = N//2+2
a = [Fraction(1, ffact(2*k)) for k in range(KK)]
s = [Fraction(0)]*KK; s[0] = Fraction(1)
for k in range(1, KK):
    s[k] = (a[k] - sum(s[i]*s[k-i] for i in range(1, k)))/2
Hf = [Fraction(0)]*(N+1)
for n in range(N+1):
    Hf[n] = 2*sum(s[k]*Fraction((-1)**(n-2*k), 2**(n-2*k)*ffact(n-2*k)) for k in range(n//2+1))
H = [mp.mpf(v.numerator)/mp.mpf(v.denominator) for v in Hf]

def Hreg(j, x):
    # H^{[j]}(x) = sum_{m>=0} H_m x^{m+j}/(m+j)!
    tot = mp.mpf(0)
    for m in range(N-j+1):
        tot += H[m]*x**(m+j)/mp.factorial(m+j)
    return tot

def kj(j, u):
    L = -mp.log(1+mp.exp(-u))
    return (1+mp.exp(-u))**(-2) * L**j / mp.factorial(j)

def BBj(j, z):
    # int_0^inf e^{-u} k_j(u) u^{-j} H^{[j]}(z u) du
    if j == 0:
        f = lambda u: mp.exp(-u)*kj(0,u)*Hreg(0, z*u)
    else:
        f = lambda u: mp.exp(-u)*kj(j,u)*u**(-j)*Hreg(j, z*u)
    return mp.quad(f, [0, 1, 4, 12, 40])

def B_exact(z, jmax=14):
    return sum(BBj(j, z) for j in range(0, jmax+1))

# ---- Pade [90/90] of B ----
P = json.load(open("pade9090.json"))
num = [mp.mpf(v) for v in P['num']]; den = [mp.mpf(v) for v in P['den']]
def peval(t):
    nv = num[0]; dv = den[0]
    for k in range(1, len(num)): nv += num[k]*t**k
    for k in range(1, len(den)): dv += den[k]*t**k
    return nv/dv

# ---- partial sums of B for |z|<rho ----
d = json.load(open("<project>/research/08_theory/D11b2_stokes_jump/step2_num450_N180.json"))
f = [mp.mpf(v) for v in d['f450']]
def B_series(z):
    tot = mp.mpf(0)
    for n in range(181):
        tot += f[n]/mp.factorial(n)*z**n
    return tot

print("=== B(zeta) on the real axis: exact Bj-sum vs Pade vs (inside radius) series ===")
for z in [mp.mpf('0.5'), mp.mpf('1.0'), mp.mpf('2.0'), mp.mpf('3.0'), mp.mpf('-1.0'), mp.mpf('-3.0'), mp.mpf('-6.0')]:
    be = B_exact(z)
    bp = peval(z)
    line = f"z={mp.nstr(z,6)}: Bexact={mp.nstr(be,25)}  Pade={mp.nstr(bp,25)}  |diff|={mp.nstr(abs(be-bp),3)}"
    if abs(z) < mp.pi/2:
        bs = B_series(z)
        line += f"  series={mp.nstr(bs,20)}  |d|={mp.nstr(abs(be-bs),3)}"
    print(line)
