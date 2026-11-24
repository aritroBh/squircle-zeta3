import sys, json, re, time
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
mp.mp.dps = 60
from hpl_series import H_series

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
R3 = json.load(open(G1 + "hpl_fit_k3_(0,0,0).json"))
def parse(name, coef):
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
    e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(","))
    num, den = (coef.split("/") + ["1"])[:2]
    return e, w, mp.mpf(int(num))/int(den)

c0 = mp.sqrt(2)/2
t0 = time.time()
F = mp.mpf(0)
for n, c in R3:
    e, w, cf = parse(n, c)
    if len(w) == 4:
        F += cf * c0**e * H_series(w, c0)
print("F =", mp.nstr(F, 30), f"({time.time()-t0:.1f}s)")
