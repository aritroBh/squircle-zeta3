
import json, re, sys, mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 800
G2 = "<project>/research/08_theory/G2_M4alpha/"
M4C = "<project>/research/08_theory/M4c_w5_closed_form/"

def hpl(w, x, N):
    a = w[-1]
    s = [mp.mpf(0)] * (N + 1)
    for m in range(1, N + 1):
        s[m] = mp.mpf(1) / m if a == 1 else mp.mpf((-1) ** (m + 1)) / m
    for a in reversed(w[:-1]):
        t = [mp.mpf(0)] * (N + 1)
        if a == 0:
            for m in range(1, N + 1): t[m] = s[m] / m
        else:
            acc = mp.mpf(0)
            for m in range(0, N):
                acc = acc * a + s[m]; t[m + 1] = acc / (m + 1)
        s = t
    return mp.fsum(s[m] * x ** m for m in range(1, N + 1))

terms5 = []
for nm, v in json.load(open(G2 + "hpl_fit_k4_(0,0,0).json")):
    m = re.match(r'c\^(-?\d+) (He?)\[(.*)\]', nm)
    if not m: continue
    w = tuple(int(t) for t in m.group(3).split(','))
    if len(w) == 5:
        terms5.append((w, int(m.group(1)), mp.mpf(Fr(v).numerator)/Fr(v).denominator))
print("weight-5 terms:", len(terms5))

c = 1/mp.sqrt(2)
N = int((760 + 25) / -mp.log10(abs(c))) + 30
print("series order N =", N)
R45 = mp.mpf(0)
cache = {}
for w, e, q in terms5:
    if w not in cache: cache[w] = hpl(w, c, N)
    nw = tuple(-a for a in w)
    if nw not in cache: cache[nw] = hpl(nw, c, N)
    R45 += q * (c**e * cache[w] + (-c)**e * cache[nw])
print("R4^{(5)}(1/sqrt2) =", mp.nstr(R45, 30))
W5 = mp.mpf(open(M4C + "W5_5400.txt").read().split()[0])
print("W5 (M4c, 5400 bits)  =", mp.nstr(W5, 30))
print("R45 - W5 =", mp.nstr(abs(R45 - W5), 3))
print("R45 + W5 =", mp.nstr(abs(R45 + W5), 3))
