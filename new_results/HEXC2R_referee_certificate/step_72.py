
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 80
import reflib as R
import json, re
from fractions import Fraction as Fr

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
psi2 = mp.mpf(0)
for key, val in json.load(open(G1 + "/hpl_fit_k3_(1,0,0).json")):
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", key)
    if not m: continue
    e = int(m.group(1))
    w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
    if len(w) != 2: continue
    psi2 += mp.mpf(str(val)) * mp.mpf(1)/2**e * R.hpl_value(w, mp.mpf(1)/2).real
L2, L3 = mp.log(2), mp.log(3)
basis = [mp.pi**2, mp.polylog(2, mp.mpf(1)/3), L2**2, L2*L3, L3**2]
names = ["pi^2", "Li2(1/3)", "log^2 2", "log2 log3", "log^2 3"]
rel = mp.pslq([psi2] + basis, tol=mp.mpf("1e-70"), maxcoeff=10**9)
print("pslq:", rel)
if rel:
    val = rel[0]*psi2 + sum(r*b for r, b in zip(rel[1:], basis))
    print("residual:", mp.nstr(abs(val), 3))
    print("=> Psi2 =", " + ".join(f"({-r}/{rel[0]})*{n}" if r % rel[0] else f"({-r//rel[0]})*{n}" for r, n in zip(rel[1:], names)))
