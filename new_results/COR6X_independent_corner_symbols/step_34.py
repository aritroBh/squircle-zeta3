t = r'''
# step_05d.py -- decisive numeric check of Omega(1/2) against the paper's closed form
import sys; sys.path.insert(0, '.')
import mpmath as mp, json, re
from hpl_num import hpl
mp.mp.dps = 150
data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[0], parse_term(s)[1], mp.mpf(c)) for s, c in data if len(parse_term(s)[1]) == 4]
def Omega(c0):
    c0 = mp.mpf(c0)
    return mp.fsum(c * (c0 ** e) * hpl(w, c0) for e, w, c in TERMS)
O = Omega(mp.mpf(1)/2)
print("Omega(1/2) =", mp.nstr(O, 50))
C = (-21*mp.polylog(4,-2) - mp.mpf(65)/8*mp.polylog(4,mp.mpf(1)/3) + mp.mpf(13)/2*mp.polylog(4,mp.mpf(2)/3)
     + mp.mpf(137)/6*mp.polylog(4,mp.mpf(1)/2) + mp.mpf(39)/64*mp.polylog(4,mp.mpf(1)/9))
print("paper combo (Li4 part) =", mp.nstr(C, 50))
D = O - C
print("D =", mp.nstr(D, 30))
l2, l3 = mp.log(2), mp.log(3)
z3 = mp.zeta(3); p2 = mp.pi**2
L2b = mp.polylog(2, mp.mpf(1)/3)
L3half = mp.polylog(3, mp.mpf(1)/2)
L3b = mp.polylog(3, mp.mpf(1)/3)
basis = [l2**4, l2**3*l3, l2**2*l3**2, l2*l3**3, l3**4,
         p2*l2**2, p2*l2*l3, p2*l3**2, mp.pi**4,
         z3*l2, z3*l3,
         L2b*l2**2, L2b*l2*l3, L2b*l3**2, L2b**2, p2*L2b,
         L3half*l2, L3b*l2, L3half*l3, L3b*l3]
vec = [D] + basis
for tol, mc in ((mp.mpf('1e-120'), 10**10),):
    rel = mp.pslq(vec, tol=tol, maxcoeff=mc, maxsteps=10000)
    print("pslq:", rel)
    if rel and rel[0] != 0:
        resid = mp.fsum(mp.mpf(r)*v for r, v in zip(rel, vec))
        print("residual recomputed:", mp.nstr(resid, 3))
        print("D/product =", -rel[0], "=> coeff list:", rel[1:])
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/step_05d.py","w").write(t)
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_05d.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=1500)
print(r.stdout); print("ERR:", (r.stderr or "")[-1000:])
