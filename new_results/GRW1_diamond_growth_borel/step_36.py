
import mpmath as mp
mp.mp.dps = 60
from fractions import Fraction
from math import factorial as ffact

N = 260
KK = N//2+2
a = [Fraction(1, ffact(2*k)) for k in range(KK)]
s_ = [Fraction(0)]*KK; s_[0] = Fraction(1)
for k in range(1, KK):
    s_[k] = (a[k] - sum(s_[i]*s_[k-i] for i in range(1, k)))/2
H = [Fraction(0)]*(N+1)
for n in range(N+1):
    H[n] = 2*sum(s_[k]*Fraction((-1)**(n-2*k), 2**(n-2*k)*ffact(n-2*k)) for k in range(n//2+1))
def Hb0(x):
    tot = mp.mpf(0); mx = mp.mpf(0)
    for m in range(N+1):
        t = mp.mpf(H[m].numerator)/mp.mpf(H[m].denominator)*x**m/mp.factorial(m)
        tot += t; mx = max(mx, abs(t))
        if m > 40 and abs(t) < mp.mpf('1e-45')*mx: break
    return tot

Hfun = lambda z: mp.sqrt(2)*mp.sqrt(1+mp.exp(-2*z))

x = mp.mpf(-20); X = -x
print("series H^[0](-20) =", mp.nstr(Hb0(x), 15))

# contour: rectangle Re s in [-1/X, X/pi+1], |Im s| <= 2X/pi+1, CCW
A = -1/X; B = X/mp.pi+1; T = 2*X/mp.pi+1
iu = mp.mpc(0,1)
def integrand(s):
    return mp.e**s * Hfun(x/s) / s
# parametrize 4 sides CCW: bottom (A-iT -> B-iT), right (B-iT -> B+iT), top (B+iT -> A+iT), left (A+iT -> A-iT)
seg1 = mp.quad(lambda sig: integrand(sig - iu*T), [A, 0, B])           # bottom
seg2 = mp.quad(lambda tau: integrand(B + iu*tau)*iu, [-T, 0, T])       # right
seg3 = mp.quad(lambda sig: integrand(sig + iu*T), [B, 0, A])           # top (B->A)
seg4 = mp.quad(lambda tau: integrand(A + iu*tau)*iu, [T, 0, -T])       # left (T->-T)
contour = (seg1+seg2+seg3+seg4)/(2*mp.pi*iu)
print("rectangle contour integral =", mp.nstr(contour, 15))
print("left side alone:", mp.nstr(seg4/(2*mp.pi*iu), 15))
print("right side alone:", mp.nstr(seg2/(2*mp.pi*iu), 15))
print("bottom:", mp.nstr(seg1/(2*mp.pi*iu), 15), " top:", mp.nstr(seg3/(2*mp.pi*iu), 15))
