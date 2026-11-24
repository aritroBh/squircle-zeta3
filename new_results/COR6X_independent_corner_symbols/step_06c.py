
# step_06c.py -- W4 integral, and identify Omega(1/sqrt2) - W4 by LLL
import sys; sys.path.insert(0, '.')
import mpmath as mp, json, re
from hpl_num import hpl
from lll_rel import find_relations
mp.mp.dps = 120
sq2 = mp.sqrt(2)
a = 1 + sq2; b = sq2 - 1
def ell1(v): return mp.log((v - b)*(v + a)/(4*v))
def ell2(v): return mp.log((a - v)*(v + b)/(4*v))
f = lambda v: ell1(v)**3 * ((8 - 7*sq2)/(48*v) - 1/(3*(v + a))) + sq2/16 * ell1(v)**2 * ell2(v) / v
W4_int = mp.quad(f, [b, a])
print("W4 integral =", mp.nstr(W4_int, 50))
print("paper W4    = 0.79578505983224081132")

# Omega(1/sqrt2)
data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[0], parse_term(s)[1], mp.mpf(c)) for s, c in data if len(parse_term(s)[1]) == 4]
def Omega(c0):
    c0 = mp.mpf(c0)
    return mp.fsum(c * (c0 ** e) * hpl(w, c0) for e, w, c in TERMS)
O = Omega(1/sq2)
print("Omega(1/sqrt2) =", mp.nstr(O, 50))
D = O - W4_int
print("Omega - W4 =", mp.nstr(D, 40))

# identify D among weight-4 constants over Q(sqrt2)
l2 = mp.log(2); L = mp.log(a); z3 = mp.zeta(3); p2 = mp.pi**2
Phi1 = -(1/(2*sq2))*mp.log(a)   # weight-1 part of Phi at 1/sqrt2: -(1-c^2)/(2c) artanh c, c=1/sqrt2
print("Phi_wt1(1/sqrt2) =", mp.nstr(Phi1, 20), " zeta3*Phi1 =", mp.nstr(z3*Phi1, 20))
names, basis = [], []
def add(nm, v): names.append(nm); basis.append(v)
add("pi4", mp.pi**4); add("sqrt2 pi4", sq2*mp.pi**4)
add("z3 l2", z3*l2); add("z3 L", z3*L); add("sqrt2 z3 l2", sq2*z3*l2); add("sqrt2 z3 L", sq2*z3*L)
add("pi2 l2^2", p2*l2**2); add("pi2 L^2", p2*L**2); add("pi2 l2 L", p2*l2*L)
add("sqrt2 pi2 l2^2", sq2*p2*l2**2); add("sqrt2 pi2 L^2", sq2*p2*L**2); add("sqrt2 pi2 l2 L", sq2*p2*l2*L)
add("l2^4", l2**4); add("L^4", L**4); add("l2^3 L", l2**3*L); add("l2 L^3", l2*L**3); add("l2^2 L^2", l2**2*L**2)
rels = find_relations([D] + basis, 100, 6)
for aa, resid in rels:
    print("   resid:", mp.nstr(resid, 3), " a0 =", aa[0], " maxcoef:", max(abs(x) for x in aa))
    if aa[0] != 0 and abs(resid) < mp.mpf('1e-80'):
        print("   REL:", [(names[i-1], c) for i, c in enumerate(aa[1:]) if c])
