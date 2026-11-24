
# step_05b.py -- numeric test: Omega(1/3) vs -(11/81)Li4(-8)-2Li4(-3)+(20/3)Li4(2/3) mod products
import sys; sys.path.insert(0, '.')
import mpmath as mp, json, re
from hpl_num import hpl
mp.mp.dps = 100

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
print("Omega(1/3) =", mp.nstr(O, 40))
print("combo      =", mp.nstr(C, 40))
print("difference =", mp.nstr(D, 40))

# product basis of weight 4 over {2,3}-periods
l2, l3 = mp.log(2), mp.log(3)
p2 = mp.pi**2
L2a, L2b, L2c = mp.polylog(2, mp.mpf(1)/2), mp.polylog(2, mp.mpf(1)/3), mp.polylog(2, mp.mpf(2)/3)
L3a, L3b, L3c = mp.polylog(3, mp.mpf(1)/2), mp.polylog(3, mp.mpf(1)/3), mp.polylog(3, mp.mpf(2)/3)
z3 = mp.zeta(3)
basis = [l2**4, l2**3*l3, l2**2*l3**2, l2*l3**3, l3**4,
         p2*l2**2, p2*l2*l3, p2*l3**2, mp.pi**4,
         z3*l2, z3*l3,
         L2a**2, L2b**2, L2c**2, L2a*L2b, L2a*L2c, L2b*L2c,
         L2a*l2**2, L2b*l2**2, L2c*l2**2, L2a*l3**2, L2b*l3**2, L2c*l3**2,
         L2a*l2*l3, L2b*l2*l3, L2c*l2*l3,
         L3a*l2, L3b*l2, L3c*l2, L3a*l3, L3b*l3, L3c*l3,
         p2*L2a, p2*L2b, p2*L2c]
vec = [D] + basis
rel = mp.pslq(vec, tol=mp.mpf('1e-90'), maxcoeff=10**6, maxsteps=10000)
print("pslq rel for [D, products]:", rel)
if rel:
    print("  D =", -mp.fsum(mp.mpf(r)*b for r, b in zip(rel[1:], basis)), " (check vs", mp.nstr(D,10), ")")
    print("  resid:", mp.nstr(rel[0]*D + mp.fsum(mp.mpf(r)*b for r, b in zip(rel[1:], basis)), 3))
