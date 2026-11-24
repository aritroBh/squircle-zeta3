
import mpmath as mp
mp.mp.dps = 90
from fractions import Fraction
from math import factorial as ffact
import json

# (a) f450 vs D11 exact zeta-forms for n=0..8
d = json.load(open("<project>/research/08_theory/D11b2_stokes_jump/step2_num450_N180.json"))
f = [mp.mpf(v) for v in d['f450']]
z2 = mp.zeta(2); z3 = mp.zeta(3); z4 = mp.zeta(4); z5 = mp.zeta(5); z6 = mp.zeta(6)
exact = {
 0: mp.mpf(1),
 1: mp.mpf(-1),
 2: z2/4 + 1,
 3: -z2/4 - z3/4 - 1,
 4: -7*z2**2/32 + z2/4 + z3/4 + z4/32 + 1,
 5: 7*z2**2/32 + 7*z2*z3/16 - z2/4 - z3/4 - z4/32 + 3*z5/16 - 1,
 6: 169*z2**3/384 - 7*z2**2/32 - 7*z2*z3/16 + 113*z2*z4/128 + z2/4 - 7*z3**2/32 + z3/4 + z4/32 - 3*z5/16 + 91*z6/192 + 1,
}
print("(a) f450 vs D11 exact forms:")
for n, v in exact.items():
    print(f"   f_{n}: residual {mp.nstr(abs(f[n]-v), 3)}")

# (b) H^[0](x) on the negative axis vs C log(2+|x|)
N = 300
KK = N//2+2
a = [Fraction(1, ffact(2*k)) for k in range(KK)]
s = [Fraction(0)]*KK; s[0] = Fraction(1)
for k in range(1, KK):
    s[k] = (a[k] - sum(s[i]*s[k-i] for i in range(1, k)))/2
H = [Fraction(0)]*(N+1)
for n in range(N+1):
    H[n] = 2*sum(s[k]*Fraction((-1)**(n-2*k), 2**(n-2*k)*ffact(n-2*k)) for k in range(n//2+1))
def Hb0(x):
    tot = mp.mpf(0); mx = mp.mpf(0)
    for m in range(N+1):
        t = mp.mpf(H[m].numerator)/mp.mpf(H[m].denominator)*x**m/mp.factorial(m)
        tot += t; mx = max(mx, abs(t))
        if m > 60 and abs(t) < mp.mpf('1e-50')*mx: break
    return tot
print("(b) |H^[0](x)| vs 40 log(2+|x|) on negative axis:")
for x in [-5, -20, -80, -200]:
    v = Hb0(mp.mpf(x))
    print(f"   H^[0]({x}) = {mp.nstr(v, 12)}   40log(2+|x|) = {mp.nstr(40*mp.log(2+abs(x)), 8)}   ok={abs(v) <= 40*mp.log(2+abs(x))}")
