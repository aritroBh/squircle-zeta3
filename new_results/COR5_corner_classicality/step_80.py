import sys, json, math, re
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
from hpl_series import hpl_coeffs

mp.mp.dps = 300
base = "<project>/research"
R3 = json.load(open(base + "/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))

def H_series2(w, x, N):
    if w == (0,): return mp.log(x)
    c = hpl_coeffs(tuple(w), N)
    tot = mp.mpf(0); xn = x
    for n in range(1, N + 1):
        tot += c[n] * xn; xn *= x
    return tot

N = int((mp.mp.dps + 15) / abs(math.log10(0.5))) + 10
c = mp.mpf(1)/2
F4_new = mp.mpf(0)
nterms = 0
for name, coef in R3:
    m = re.match(r'c\^(-?\d+) H\[(.*)\]', name)
    e = int(m.group(1))
    w = tuple(int(z) for z in m.group(2).split(','))
    if len(w) != 4: continue
    num, den = (coef.split("/") + ["1"])[:2]
    F4_new += mp.mpf(int(num))/int(den) * c**e * H_series2(w, c, N)
    nterms += 1
F4_old = mp.mpf(open("hex_F4.txt").read())
print("terms:", nterms)
print("F4 fresh:", mp.nstr(F4_new, 40))
print("F4 saved:", mp.nstr(F4_old, 40))
print("|diff|:", mp.nstr(abs(F4_new - F4_old), 3))
