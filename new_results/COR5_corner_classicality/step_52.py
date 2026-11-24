import sys, json, math, time
sys.path.insert(0, "<project>/research/10_frontier/G1_zeta3_any_angle")
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
from hpl_series import H_series
import re

base = "<project>/research"
R3 = json.load(open(base + "/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
Psi = json.load(open(base + "/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(1,0,0).json"))
def parse(name, coef):
    if name == '1': return 0, (), mp.mpf(1)
    if name == 'c': return 1, (), mp.mpf(1)
    m = re.match(r'c\^(-?\d+) H\[(.*)\]', name)
    num, den = (coef.split("/") + ["1"])[:2]
    return int(m.group(1)), tuple(int(z) for z in m.group(2).split(',')), mp.mpf(int(num))/int(den)

def combo(terms, cc, N):
    tot = mp.mpf(0)
    for name, coef in terms:
        e, w, cf = parse(name, coef)
        val = mp.mpf(1) if w == () else H_series(w, cc, N)
        tot += cf * cc**e * val
    return tot

def Phi(c):
    return mp.mpf(1)/2 - (1 - c*c)/(2*c)*mp.atanh(c)

mp.mp.dps = 760
N = int((mp.mp.dps + 15) / abs(math.log10(0.5))) + 10
print("series length N =", N)
t0 = time.time()
c = mp.mpf(1)/2
val = -mp.zeta(3)*Phi(c) + mp.zeta(2)*combo(Psi, c, N) + combo(R3, c, N)
t1 = time.time()
print(f"HPL formula at dps 760: {t1-t0:.1f}s")
v1100 = mp.mpf(open("M3_2pi3_1100.txt").read())
print("M3 HPL formula:", mp.nstr(val, 60))
print("|HPL - quadrature1100| =", mp.nstr(abs(val - v1100), 3))
open("M3_2pi3_hpl_760.txt", "w").write(mp.nstr(val, 755))
