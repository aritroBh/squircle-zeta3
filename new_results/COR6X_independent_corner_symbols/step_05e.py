
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
O2 = Omega(mp.mpf(1)/2)
C2 = (-21*mp.polylog(4,-2) - mp.mpf(65)/8*mp.polylog(4,mp.mpf(1)/3) + mp.mpf(13)/2*mp.polylog(4,mp.mpf(2)/3)
     + mp.mpf(137)/6*mp.polylog(4,mp.mpf(1)/2) + mp.mpf(39)/64*mp.polylog(4,mp.mpf(1)/9))
D2 = O2 - C2
print("D2 =", mp.nstr(D2, 40))
l2, l3 = mp.log(2), mp.log(3)
z3 = mp.zeta(3); p2 = mp.pi**2
L2b = mp.polylog(2, mp.mpf(1)/3)
L3half = mp.polylog(3, mp.mpf(1)/2); L3b = mp.polylog(3, mp.mpf(1)/3)
basis = [l2**4, l2**3*l3, l2**2*l3**2, l2*l3**3, l3**4,
         p2*l2**2, p2*l2*l3, p2*l3**2, mp.pi**4,
         z3*l2, z3*l3,
         L2b*l2**2, L2b*l2*l3, L2b*l3**2, L2b**2, p2*L2b,
         L3half*l2, L3b*l2, L3half*l3, L3b*l3]
print("LLL relations for [D2, products] at 180 digits:")
for a, resid, last in find_relations([D2] + basis, 180, 8):
    print("   vec[:4]... resid:", mp.nstr(resid, 3), " a0 =", a[0], " small coeffs:", max(abs(c) for c in a) < 1000)
    if resid == 0 or abs(resid) < mp.mpf('1e-150'):
        print("     RELATION:", a)
