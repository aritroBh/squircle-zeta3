# Planner check: the number Omega(1/2) (the 56 weight-4 HPL terms of R3, G1's JSON) against COR5's closed form F4 (5 Li4 + products).
import json, re, mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 140
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
def hpl(w, x, N):
    a = w[-1]; s = [mp.mpf(0)] * (N + 1)
    for m in range(1, N + 1): s[m] = mp.mpf(1) / m if a == 1 else mp.mpf((-1) ** (m + 1)) / m
    for a in reversed(w[:-1]):
        t = [mp.mpf(0)] * (N + 1)
        if a == 0:
            for m in range(1, N + 1): t[m] = s[m] / m
        else:
            acc = mp.mpf(0)
            for m in range(0, N): acc = acc * a + s[m]; t[m + 1] = acc / (m + 1)
        s = t
    return mp.fsum(s[m] * x ** m for m in range(1, N + 1))
c = mp.mpf(1) / 2; N = 520; n4 = 0; Om = mp.mpf(0)
for nm, v in json.load(open(G1 + "hpl_fit_k3_(0,0,0).json")):
    if nm == '1': continue
    m = re.match(r'c\^(-?\d+) (He?)\[(.*)\]', nm); w = tuple(int(t) for t in m.group(3).split(','))
    if len(w) != 4: continue
    n4 += 1; q = mp.mpf(Fr(v).numerator) / Fr(v).denominator; e = int(m.group(1))
    f = c ** e * hpl(w, c, N)
    if m.group(2) == 'He': f += (-c) ** e * hpl(w, -c, N)
    Om += q * f
Li = lambda n, x: mp.polylog(n, x); l2, l3 = mp.log(2), mp.log(3); pi2 = mp.pi ** 2; z3 = mp.zeta(3); Q = mp.mpf; L213 = Li(2, Q(1)/3)
F4 = (-21*Li(4, -2) - Q(65)/8*Li(4, Q(1)/3) + Q(13)/2*Li(4, Q(2)/3) + Q(137)/6*Li(4, Q(1)/2) + Q(39)/64*Li(4, Q(1)/9)
      + Q(15)/2*Li(3, Q(1)/3)*l2 + Q(37)/2*Li(3, -Q(1)/2)*l2 - Q(5)/4*Li(3, -Q(1)/2)*l3 + Q(119)/16*z3*l2 - Q(37)/32*z3*l3
      + L213**2 - Q(23)/4*L213*l2**2 + Q(7)/4*L213*l2*l3 + L213*l3**2 - pi2/6*L213
      - Q(455)/144*pi2*l2**2 + Q(9)/8*pi2*l2*l3 - Q(7)/32*pi2*l3**2 - Q(949)/2160*mp.pi**4
      - Q(211)/144*l2**4 + Q(101)/24*l2**3*l3 - Q(15)/8*l2**2*l3**2 - Q(35)/24*l2*l3**3 + Q(113)/192*l3**4)
li4part = -21*Li(4, -2) - Q(65)/8*Li(4, Q(1)/3) + Q(13)/2*Li(4, Q(2)/3) + Q(137)/6*Li(4, Q(1)/2) + Q(39)/64*Li(4, Q(1)/9)
print("weight-4 terms used:", n4)
print("Omega(1/2) =", mp.nstr(Om, 40))
print("Omega(1/2) - five-Li4 part =", mp.nstr(Om - li4part, 20), " (COR6X quotes -51.777095579...)")
print("|Omega(1/2) - COR5 closed form F4| =", mp.nstr(abs(Om - F4), 3), " at dps", mp.mp.dps)
