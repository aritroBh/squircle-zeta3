
# step_05g.py -- comprehensive product basis, decisive for Omega(1/2) and Omega(1/3)
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
Xs = [mp.mpf(1)/2, mp.mpf(1)/3, mp.mpf(2)/3, mp.mpf(1)/4, mp.mpf(3)/4, mp.mpf(1)/6,
      mp.mpf(1)/8, mp.mpf(1)/9, mp.mpf(2)/9, -mp.mpf(1)/2, -mp.mpf(1)/3, -mp.mpf(2)/3]
L2 = {str(x): mp.polylog(2, x) for x in Xs}
L3 = {str(x): mp.polylog(3, x) for x in Xs}
names, basis = [], []
def add(name, val):
    names.append(name); basis.append(val)
add("l2^4", l2**4); add("l2^3 l3", l2**3*l3); add("l2^2 l3^2", l2**2*l3**2)
add("l2 l3^3", l2*l3**3); add("l3^4", l3**4)
add("pi2 l2^2", p2*l2**2); add("pi2 l2 l3", p2*l2*l3); add("pi2 l3^2", p2*l3**2); add("pi4", mp.pi**4)
add("z3 l2", z3*l2); add("z3 l3", z3*l3)
for x in Xs:
    sx = str(x)
    add(f"pi2 Li2({sx})", p2*L2[sx])
for x in Xs:
    sx = str(x)
    add(f"Li2({sx}) l2^2", L2[sx]*l2**2); add(f"Li2({sx}) l2 l3", L2[sx]*l2*l3); add(f"Li2({sx}) l3^2", L2[sx]*l3**2)
for i, x in enumerate(Xs):
    for y in Xs[i:]:
        add(f"Li2({x})Li2({y})", L2[str(x)]*L2[str(y)])
for x in Xs:
    sx = str(x)
    add(f"Li3({sx}) l2", L3[sx]*l2); add(f"Li3({sx}) l3", L3[sx]*l3)
print("basis size:", len(basis))

print("================ c0 = 1/2")
O2 = Omega(mp.mpf(1)/2)
C2 = (-21*mp.polylog(4,-2) - mp.mpf(65)/8*mp.polylog(4,mp.mpf(1)/3) + mp.mpf(13)/2*mp.polylog(4,mp.mpf(2)/3)
     + mp.mpf(137)/6*mp.polylog(4,mp.mpf(1)/2) + mp.mpf(39)/64*mp.polylog(4,mp.mpf(1)/9))
D2 = O2 - C2
print("D2 =", mp.nstr(D2, 25))
rels = find_relations([D2] + basis, 180, 10)
found = False
for a, resid in rels:
    if a[0] != 0 and abs(resid) < mp.mpf('1e-140'):
        found = True
        print("  SMALL RELATION involving D2:", a, " resid:", mp.nstr(resid,3))
        print("   terms:", [(names[i-1], c) for i, c in enumerate(a[1:]) if c])
if not found:
    print("  NO relation involving D2 with small coefficients (checked 10 shortest lattice vectors):")
    for a, resid in rels[:6]:
        print("    a0 =", a[0], " maxcoef ~", max(abs(c) for c in a), " resid:", mp.nstr(resid, 2))

print("================ c0 = 1/3")
O3 = Omega(mp.mpf(1)/3)
C3 = (-mp.mpf(11)/81 * mp.polylog(4, -8) - 2 * mp.polylog(4, -3) + mp.mpf(20)/3 * mp.polylog(4, mp.mpf(2)/3))
D3 = O3 - C3
print("D3 =", mp.nstr(D3, 25))
rels = find_relations([D3] + basis, 180, 10)
found = False
for a, resid in rels:
    if a[0] != 0 and abs(resid) < mp.mpf('1e-140'):
        found = True
        print("  SMALL RELATION involving D3:", a, " resid:", mp.nstr(resid,3))
        print("   terms:", [(names[i-1], c) for i, c in enumerate(a[1:]) if c])
if not found:
    print("  NO relation involving D3 with small coefficients:")
    for a, resid in rels[:6]:
        print("    a0 =", a[0], " maxcoef ~", max(abs(c) for c in a), " resid:", mp.nstr(resid, 2))
