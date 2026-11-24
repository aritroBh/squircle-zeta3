
import sys; sys.path.insert(0, '.')
import mpmath as mp, json, re
from hpl_num import hpl
mp.mp.dps = 100
l2, l3 = mp.log(2), mp.log(3)
L2a, L2b, L2c = mp.polylog(2, mp.mpf(1)/2), mp.polylog(2, mp.mpf(1)/3), mp.polylog(2, mp.mpf(2)/3)
# check the exact product relations used for reduction
print("Li2(1/2) = pi^2/12 - l2^2/2:", mp.nstr(L2a - (mp.pi**2/12 - l2**2/2), 3))
print("Li2(2/3) = pi^2/6 + l2*l3 - l3^2 - Li2(1/3):", mp.nstr(L2c - (mp.pi**2/6 + l2*l3 - l3**2 - L2b), 3))
# the pslq relation from before:
chk = L2b**2 + L2a*L2b - 2*L2a*L2c - L2a*l2**2 + L2a*l3**2 - L2a*l2*l3
print("previous pslq relation residue:", mp.nstr(chk, 3))

data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[0], parse_term(s)[1], mp.mpf(c)) for s, c in data if len(parse_term(s)[1]) == 4]
def Omega(c0):
    c0 = mp.mpf(c0)
    return mp.fsum(c * (c0 ** e) * hpl(w, c0) for e, w, c in TERMS)
O = Omega(mp.mpf(1)/3)
C = (-mp.mpf(11)/81 * mp.polylog(4, -8) - 2 * mp.polylog(4, -3) + mp.mpf(20)/3 * mp.polylog(4, mp.mpf(2)/3))
D = O - C
print("D =", mp.nstr(D, 30))
z3 = mp.zeta(3); p2 = mp.pi**2
basis = [l2**4, l2**3*l3, l2**2*l3**2, l2*l3**3, l3**4,
         p2*l2**2, p2*l2*l3, p2*l3**2, mp.pi**4,
         z3*l2, z3*l3,
         L2b*l2**2, L2b*l2*l3, L2b*l3**2, L2b**2, p2*L2b]
vec = [D] + basis
rel = mp.pslq(vec, tol=mp.mpf('1e-85'), maxcoeff=10**8, maxsteps=10000)
print("pslq on reduced basis:", rel)
if rel and rel[0] != 0:
    print("D as product:", [r for r in rel])
