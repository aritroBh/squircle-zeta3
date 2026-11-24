
import mpmath as mp
import json, time
mp.mp.dps = 400

# exact H_n
from fractions import Fraction
from math import factorial as ffact
N = 182
KK = N//2+2
a = [Fraction(1, ffact(2*k)) for k in range(KK)]
s = [Fraction(0)]*KK; s[0] = Fraction(1)
for k in range(1, KK):
    s[k] = (a[k] - sum(s[i]*s[k-i] for i in range(1, k)))/2
Hf = [Fraction(0)]*(N+1)
for n in range(N+1):
    Hf[n] = 2*sum(s[k]*Fraction((-1)**(n-2*k), 2**(n-2*k)*ffact(n-2*k)) for k in range(n//2+1))
H = [mp.mpf(v.numerator)/mp.mpf(v.denominator) for v in Hf]

d = json.load(open("<project>/research/08_theory/D11b2_stokes_jump/step2_num450_N180.json"))
f = [mp.mpf(v) for v in d['f450']]
b = [f[n]/mp.factorial(n) for n in range(len(f))]
rser = [b[n]-H[n] for n in range(181)]
print("radius check |r_n|^(1/n) vs 1/(2 rho)=1/pi:")
for n in [90, 140, 180]:
    print(f"  n={n}: {mp.nstr(abs(rser[n])**(mp.mpf(1)/n), 12)}   1/pi={mp.nstr(1/mp.pi,12)}")

def pade_coeffs(a, L, M):
    A = mp.matrix(M, M); rhs = mp.matrix(M, 1)
    for i in range(M):
        for j in range(M):
            k = L + i - j
            A[i, j] = a[k] if k >= 0 else mp.mpf(0)
        rhs[i] = -a[L + 1 + i]
    q = (A**-1)*rhs
    den = [mp.mpf(1)] + [q[j] for j in range(M)]
    num = [sum(den[j]*a[n-j] for j in range(0, min(M, n)+1)) for n in range(L+1)]
    return num, den
t0=time.time()
num, den = pade_coeffs(rser, 90, 90)
print("R-pade build s:", round(time.time()-t0,1))
# certification
ser = [mp.mpf(0)]*181
ser[0] = num[0]
for n in range(1, 181):
    ser[n] = (num[n] if n < len(num) else mp.mpf(0)) - sum(den[j]*ser[n-j] for j in range(1, min(n, 90)+1))
worst = max(abs(ser[n]-rser[n]) for n in range(181))
print("R-pade reproduces r_0..r_180, max err:", mp.nstr(worst, 4))
json.dump({"num":[mp.nstr(v,380) for v in num], "den":[mp.nstr(v,380) for v in den]}, open("padeR9090.json","w"))
