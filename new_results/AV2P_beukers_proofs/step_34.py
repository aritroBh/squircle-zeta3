from mpmath import mp, mpf, zeta, log, pi, nprint, quad, clsin, catalan
mp.dps = 60
G = catalan; Cl = clsin(2, pi/3)
odd = mpf(open('odd2.txt').read())
# J* = int_0^1 [arcsin(sqrt(t)/2)/t - arcsin(t/2)]/(1-t) dt ; odd = (pi/4) J*  ??
Js = quad(lambda t: (mp.asin(mp.sqrt(t)/2)/t - mp.asin(t/2))/(1-t), [0, 1])
print("J* =", nprint(Js, 40))
print("(pi/4)J* =", nprint(pi/4*Js, 40))
print("odd =", nprint(odd, 40))
print("odd - (pi/4)J*:", nprint(odd - pi/4*Js, 5))
print("J* - 4G/3:", nprint(Js - 4*G/3, 5))
print("J* - (4G/3 - pi log2/6):", nprint(Js - (4*G/3 - pi*log(2)/6), 5))
# analytic pieces
first = 2*quad(lambda th: th/tan(th), [0, pi/6]) if False else None
F1 = 2*quad(lambda th: th*mp.cos(th)/mp.sin(th), [0, pi/6])
print("2 int th cot =", nprint(F1, 30), " Cl =", nprint(Cl, 30))
F2 = 2*quad(lambda th: th*mp.cos(th)/(1+2*mp.sin(th)), [0, pi/6])
print("2 int th cos/(1+2sin) =", nprint(F2, 30))
print("F1 - F2 =", nprint(F1-F2, 40))
