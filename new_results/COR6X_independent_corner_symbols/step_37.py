t = r'''
# step_05f.py -- decisive numeric tests of the paper's Omega closed forms vs products
import sys; sys.path.insert(0, '.')
import mpmath as mp, json, re
from hpl_num import hpl
from lll_rel import find_relations
mp.mp.dps = 200
data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[0], parse_term(s)[1], mp.mpf(c)) for s, c in data if len(parse_term(s)[1]) == 4]
def Omega(c0):
    c0 = mp.mpf(c0)
    return mp.fsum(c * (c0 ** e) * hpl(w, c0) for e, w, c in TERMS)

l2, l3 = mp.log(2), mp.log(3)
z3 = mp.zeta(3); p2 = mp.pi**2
L2b = mp.polylog(2, mp.mpf(1)/3)
L3half = mp.polylog(3, mp.mpf(1)/2); L3b = mp.polylog(3, mp.mpf(1)/3)
basis = [l2**4, l2**3*l3, l2**2*l3**2, l2*l3**3, l3**4,
         p2*l2**2, p2*l2*l3, p2*l3**2, mp.pi**4,
         z3*l2, z3*l3,
         L2b*l2**2, L2b*l2*l3, L2b*l3**2, L2b**2, p2*L2b,
         L3half*l2, L3b*l2, L3half*l3, L3b*l3]

print("================ c0 = 1/2")
O2 = Omega(mp.mpf(1)/2)
C2 = (-21*mp.polylog(4,-2) - mp.mpf(65)/8*mp.polylog(4,mp.mpf(1)/3) + mp.mpf(13)/2*mp.polylog(4,mp.mpf(2)/3)
     + mp.mpf(137)/6*mp.polylog(4,mp.mpf(1)/2) + mp.mpf(39)/64*mp.polylog(4,mp.mpf(1)/9))
D2 = O2 - C2
print("D2 = Omega(1/2) - paper Li4 combo =", mp.nstr(D2, 30))
rels = find_relations([D2] + basis, 180, 6)
for a, resid in rels:
    tag = "RELATION" if abs(resid) < mp.mpf('1e-150') else ""
    print("   resid:", mp.nstr(resid, 3), " a0 =", a[0], " maxcoef =", max(abs(c) for c in a), tag)
    if tag: print("     ", a)

print("================ c0 = 1/3")
O3 = Omega(mp.mpf(1)/3)
C3 = (-mp.mpf(11)/81 * mp.polylog(4, -8) - 2 * mp.polylog(4, -3) + mp.mpf(20)/3 * mp.polylog(4, mp.mpf(2)/3))
D3 = O3 - C3
print("D3 = Omega(1/3) - paper Li4 combo =", mp.nstr(D3, 30))
rels = find_relations([D3] + basis, 180, 6)
for a, resid in rels:
    tag = "RELATION" if abs(resid) < mp.mpf('1e-150') else ""
    print("   resid:", mp.nstr(resid, 3), " a0 =", a[0], " maxcoef =", max(abs(c) for c in a), tag)
    if tag: print("     ", a)
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/step_05f.py","w").write(t)
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_05f.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=1500)
print(r.stdout); print("ERR:", (r.stderr or "")[-1000:])
