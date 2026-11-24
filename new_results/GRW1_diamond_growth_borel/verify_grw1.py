#!/usr/bin/env python
# verify_grw1.py - consolidated verification for GRW1 (target: < 5 min)
# Labels: [EXACT] exact rational arithmetic, [NUM] mpmath numerical, [DATA] D11b2 450-digit table, [CERT] arb
import json, time
import mpmath as mp
from fractions import Fraction
from math import factorial as ffact
mp.mp.dps = 60
T0 = time.time()
def elapsed(): return f"{time.time()-T0:6.1f}s"
print("="*78)
print("GRW1 verification: diamond-end growth law and Borel summability", elapsed())
print("="*78)

# ---------------------------------------------------------------- [EXACT] H_n
print("\n[1] Exact rational H_n (square), H(z)=2e^{-z/2}sqrt(cosh z)=sqrt(2(1+e^{-2z}))")
N = 190
KK = N//2+2
a = [Fraction(1, ffact(2*k)) for k in range(KK)]
s = [Fraction(0)]*KK; s[0] = Fraction(1)
for k in range(1, KK):
    s[k] = (a[k] - sum(s[i]*s[k-i] for i in range(1, k)))/2
Hf = [Fraction(0)]*(N+1)
for n in range(N+1):
    Hf[n] = 2*sum(s[k]*Fraction((-1)**(n-2*k), 2**(n-2*k)*ffact(n-2*k)) for k in range(n//2+1))
ok = True
for n in range(0, N+1):
    lhs = sum(Hf[k]*Hf[n-k] for k in range(n+1))
    rhs = Fraction(4) if n == 0 else 2*Fraction((-2)**n, ffact(n))
    if lhs != rhs: ok = False; break
print("    H^2 == 2(1+e^{-2z}) as formal series to order 190:", ok, "[EXACT]")
H = [mp.mpf(v.numerator)/mp.mpf(v.denominator) for v in Hf]

# ------------------------------------------------- [EXACT+NUM] Darboux laws
print("\n[2] Growth law for H_n: lawH(n) = sqrt2 cos(pi(3-2n)/4) (2/pi)^n n^-3/2")
twopi = 2/mp.pi
def lawH(n):
    return mp.sqrt(2)*mp.cos(mp.pi*(3-2*n)/4)*twopi**n*n**(-mp.mpf('1.5'))
worst = 0; C2 = 0; minr = 99
for n in range(1, N+1):
    r = abs(H[n]/lawH(n))
    worst = max(worst, r)
    law2 = lawH(n)*(1 + (3 - 3*mp.pi*((-1)**n))/(8*n))
    C2 = max(C2, abs(H[n]-law2)/(twopi**n*n**(-mp.mpf('2.5'))))
    if n >= 10: minr = min(minr, abs(H[n]/lawH(n)))
print(f"    max |H_n| (pi/2)^n n^1.5 for n<=190 = {mp.nstr(worst, 12)}  (explicit C_H in |H_n| <= C_H (2/pi)^n n^-3/2)")
print(f"    min |H_n/lawH(n)| for 10<=n<=190 = {mp.nstr(minr, 12)}  (explicit lower constant)")
print(f"    |H_n - lawH(1+(3-3pi(-1)^n)/(8n))| <= {mp.nstr(C2,6)} (2/pi)^n n^-5/2  for n<=190")
for n in [40, 80, 120, 180]:
    print(f"    n={n}: H_n/lawH = {mp.nstr(H[n]/lawH(n), 22)}")

# ------------------------------------------------- [DATA] f_n
print("\n[3] f_n (D11b2 450-digit table) vs growth law and vs n! H_n")
d = json.load(open("<project>/research/08_theory/D11b2_stokes_jump/step2_num450_N180.json"))
f = [mp.mpf(v) for v in d['f450']]
def law(n):
    return lawH(n)*mp.factorial(n)
for n in [40, 80, 120, 160, 180]:
    print(f"    n={n}: f_n/law = {mp.nstr(f[n]/law(n), 22)}   f_n/(n! H_n) = {mp.nstr(f[n]/(mp.factorial(n)*H[n]), 22)}")
import numpy as np
for r in range(4):
    ns = [n for n in range(100, 181) if n % 4 == r]
    X = np.array([[1.0, 1.0/n, 1.0/n**2] for n in ns])
    Y = np.array([float(f[n]/law(n)) for n in ns])
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    print(f"    class {r} mod 4: Richardson const = {coef[0]:.10f}")

# ------------------------------------------------- [EXACT] tail bound
print("\n[4] Uniformity step: |f_n - n! H_n eta(n)| <= tail, with explicit bound")
CH = worst
rho = mp.pi/2
for n in [20, 40, 80, 120]:
    tail_exact = mp.mpf(0)
    for j in range(1, n+1):
        tail_exact += abs(H[n-j])*mp.factorial(n-j)/(mp.factorial(j)*(j+1)**(n-j+1))
    eta_n = (1-mp.mpf(2)**(1-n))*mp.zeta(n)
    main = abs(H[n])*mp.factorial(n)*eta_n
    bound = CH*mp.e**(2*rho)*rho**(-n)*mp.factorial(n)*mp.mpf(2)**(-n-1)/n
    print(f"    n={n}: tail/main={mp.nstr(tail_exact/main, 4)}   explicit bound/main={mp.nstr(bound/main, 4)}")

# ------------------------------------------------- [NUM] identity
print("\n[5] Identity I(eta) = int_0^inf e^{-u} H(eta u) (1+e^{-u})^{-2-eta} du vs direct")
Hfun = lambda z: mp.sqrt(2)*mp.sqrt(1+mp.exp(-2*z))
for eta in ['0.3','-0.3']:
    e = mp.mpf(eta)
    Iu = mp.quad(lambda u: mp.exp(-u)*Hfun(e*u)*(1+mp.exp(-u))**(-2-e), [0, 1, 5, 20, 100])
    It = mp.quad(lambda t: (t*(1-t))**(e/2)*mp.sqrt(mp.cosh(e*mp.log(t/(1-t)))), [0, mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('0.75'), 1])
    print(f"    eta={eta}: |Iu - It| = {mp.nstr(abs(Iu-It), 3)}")

# ------------------------------------------------- [NUM+CERT] Borel spot check
print("\n[6] Borel summability: Pade-Borel [30/30] from f_0..f_60 at eta=0.1, vs CERTIFIED reference")
old = mp.mp.dps; mp.mp.dps = 220
f220 = [mp.mpf(v) for v in d['f450']]
b = [f220[n]/mp.factorial(n) for n in range(61)]
M = 30; L = 30
A = mp.matrix(M, M); rhs = mp.matrix(M, 1)
for i in range(M):
    for j in range(M):
        A[i, j] = b[L+i-j]
    rhs[i] = -b[L+1+i]
q = (A**-1)*rhs
den = [mp.mpf(1)] + [q[j] for j in range(M)]
num = [sum(den[j]*b[n-j] for j in range(0, min(M, n)+1)) for n in range(L+1)]
def peval(t):
    nv = num[0]; dv = den[0]
    for k in range(1, len(num)): nv += num[k]*t**k
    for k in range(1, len(den)): dv += den[k]*t**k
    return nv/dv
mp.mp.dps = old
S = mp.quad(lambda t: mp.exp(-t)*peval(mp.mpf('0.1')*t), [0, 2, 8, 20, 50])
cert = json.load(open("iref_certified.json"))
ref = mp.mpf(cert['0.1'][0].split(' +/-')[0].strip('[]'))
print(f"    S(0.1)={mp.nstr(S, 40)}")
print(f"    |S - certified I(0.1)| = {mp.nstr(abs(S-ref), 3)}  ({mp.nstr(-mp.log10(abs(S-ref)),4)} digits)")

# ------------------------------------------------- singular directions
print("\n[7] Singular (Stokes) directions of the Borel transform")
print("    square (L=0): arg zeta = +/- pi/2 (tower z = i pi (k+1/2)); real axis clear")
import math
Lval = math.log(2)
ths = [math.atan2(math.pi*(2*k+1)/2, Lval) for k in range(4)]
print(f"    rhombus L=log2: theta_k = {[round(t,6) for t in ths]} (+ pi - theta_k); 0 and pi are NOT Stokes")

# ------------------------------------------------- [NUM] Delta constants E*
print("\n[8] Delta-domain constant E* = max |E(is)| on the slit (should be <= 5/24)")
def Emp(s):
    w = 1j*s
    return abs((mp.sqrt((1-mp.e**(-2*w))/(2*w)) - 1 + w/2)/(w*w))
mx = max(Emp(i*mp.pi/2/1000) for i in range(1, 1001))
print(f"    max |E(is)| on (0,pi/2] (1000 pts) = {mp.nstr(mx, 15)}   5/24 = {mp.nstr(mp.mpf(5)/24, 15)}  ok={mx <= mp.mpf(5)/24}")

# -------------------------------------- [EXACT] rhombus H_n (q=5/13) vs law
print("\n[9] Rhombus exact H_n (q = 5/13, A=3,B=2) vs the two-term Darboux law")
qR = Fraction(5, 13)
NR = 120
aR = [(Fraction(1, ffact(n)) if n%2==0 else qR/Fraction(ffact(n))) for n in range(NR+1)]
gR = [Fraction(0)]*(NR+1); gR[0] = Fraction(1)
for k in range(1, NR+1):
    gR[k] = (aR[k] - sum(gR[i]*gR[k-i] for i in range(1, k)))/2
HRf = [Fraction(0)]*(NR+1)
for n in range(NR+1):
    HRf[n] = sum((gR[k] + ((-1)**k)*gR[k]) * Fraction((-1)**(n-k), 2**(n-k)*ffact(n-k)) for k in range(0, n+1))
HR = [mp.mpf(v.numerator)/mp.mpf(v.denominator) for v in HRf]
A_, B_ = 3, 2
r_ = mp.sqrt(13); L_ = mp.log(mp.mpf(2)/3); rho_ = mp.sqrt(L_**2 + mp.pi**2/4)
th0_ = mp.atan2(mp.pi/2, L_)
zl = [L_ + 1j*mp.pi/2, -L_ + 1j*mp.pi/2, L_ - 1j*mp.pi/2, -L_ - 1j*mp.pi/2]
cl = [A_*mp.sqrt(2)/r_, B_*mp.sqrt(2)/r_, A_*mp.sqrt(2)/r_, B_*mp.sqrt(2)/r_]
def bracket_(n):
    return B_*((-1)**n)*mp.cos((n-mp.mpf('0.5'))*th0_) - A_*mp.sin((n-mp.mpf('0.5'))*th0_)
def law2R(n):
    tot = mp.mpc(0)
    for z, c in zip(zl, cl):
        tot += c*mp.sqrt(-z)*z**(-n)*(1 + mp.mpf(3)/(8*n) - 3*z/(4*n))
    return (-1/(2*mp.sqrt(mp.pi)))*n**(-mp.mpf('1.5'))*tot
errs = [abs(1-HR[n]/law2R(n).real) for n in range(60, NR+1) if abs(bracket_(n)) > 0.6]
print(f"    worst |1 - H_n/law2H| at safe n in [60,{NR}]: {mp.nstr(max(errs), 4)}")
print("\nDONE", elapsed())

