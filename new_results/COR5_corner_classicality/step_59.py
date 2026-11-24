import sys, json, math, time
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
from hpl_series import H_series
import re

base = "<project>/research"
R3 = json.load(open(base + "/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
Psi = json.load(open(base + "/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(1,0,0).json"))
def parse(name, coef):
    num, den = (coef.split("/") + ["1"])[:2]
    cf = mp.mpf(int(num))/int(den)
    if name == '1': return 0, (), cf
    if name == 'c': return 1, (), cf
    m = re.match(r'c\^(-?\d+) H\[(.*)\]', name)
    return int(m.group(1)), tuple(int(z) for z in m.group(2).split(',')), cf

def combo(terms, cc, N, weights=None):
    tot = mp.mpf(0)
    for name, coef in terms:
        e, w, cf = parse(name, coef)
        if weights is not None and len(w) not in weights: continue
        val = mp.mpf(1) if w == () else H_series(w, cc, N)
        tot += cf * cc**e * val
    return tot

def Phi(c):
    return mp.mpf(1)/2 - (1 - c*c)/(2*c)*mp.atanh(c)

mp.mp.dps = 1120
N = int((mp.mp.dps + 15) / abs(math.log10(0.5))) + 10
c = mp.mpf(1)/2
t0 = time.time()
# full M3
M3full = -mp.zeta(3)*Phi(c) + mp.zeta(2)*combo(Psi, c, N) + combo(R3, c, N)
# weight-4 part F = R3^{(4)}(1/2)
F4 = combo(R3, c, N, weights={4})
# weight-3 part of R3
R3_3 = combo(R3, c, N, weights={3})
# weight-2 part of R3
R3_2 = combo(R3, c, N, weights={2})
# Psi parts
Psi2 = combo(Psi, c, N, weights={2})
Psi1 = combo(Psi, c, N, weights={1})
Psi0 = combo(Psi, c, N, weights={0})
t1 = time.time()
print(f"computed in {t1-t0:.1f}s, N={N}")
v1100 = mp.mpf(open("M3_2pi3_1100.txt").read())
print("|M3 HPL(1120) - quad(1100)| =", mp.nstr(abs(M3full - v1100), 3))
open("M3_2pi3_hpl_1120.txt", "w").write(mp.nstr(M3full, 1115))
# weight decomposition: M3 = W4 + W3 + W2
atanh = mp.atanh(c)
W4 = F4 + mp.zeta(2)*Psi2 + mp.zeta(3)*((1-c*c)/(2*c))*atanh
W3 = R3_3 + mp.zeta(2)*Psi1 - mp.zeta(3)*mp.mpf(1)/2
W2 = R3_2 + mp.zeta(2)*mp.mpf(1)/2
print("W4+W3+W2 == M3:", mp.nstr(abs(W4+W3+W2 - M3full), 3))
print("Phi check: -z3*Phi = -z3/2 + z3*((1-c^2)/(2c))*atanh:",
      mp.nstr(abs(-mp.zeta(3)*Phi(c) - (-mp.zeta(3)/2 + mp.zeta(3)*((1-c*c)/(2*c))*atanh)), 3))
print("Psi0 check:", mp.nstr(Psi0, 10), "(expect 0.5)")
for nm, v in [("F4", F4), ("R3_3", R3_3), ("R3_2", R3_2), ("Psi1", Psi1), ("Psi2", Psi2),
              ("W4", W4), ("W3", W3), ("W2", W2)]:
    open(f"hex_{nm}.txt", "w").write(mp.nstr(v, 1115))
    print(nm, "=", mp.nstr(v, 40))
