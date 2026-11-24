"""step_29.py -- derive the ladder H(z;xi) = int_0^z Li2(w)/(w-xi) dw numerically (small LLL
for identity derivation), then verify by numerical differentiation at random points.
"""
import mpmath as mp
from ir_utils import find_relations
mp.mp.dps = 100

def Hval(z, xi):
    return mp.quad(lambda w: mp.polylog(2, w)/(w - xi), [0, z/2, z])

xi = mp.mpf('2.7'); z = mp.mpf('0.31')
Hz = Hval(z, xi)
L = mp.log
cand = []
def add(name, val):
    cand.append((name, val))
add("Li3(z)", mp.polylog(3, z))
add("Li3(z/xi)", mp.polylog(3, z/xi))
add("Li3(1-z)", mp.polylog(3, 1-z))
add("Li3((xi-z)/xi)", mp.polylog(3, (xi-z)/xi))
add("Li3(z/(z-xi))", mp.polylog(3, z/(z-xi)))
add("Li3(xi/(xi-z))", mp.polylog(3, xi/(xi-z)))
add("Li3((z-xi)/(z-1))", mp.polylog(3, (z-xi)/(z-1)))
add("Li3((xi-z)/(xi-1))", mp.polylog(3, (xi-z)/(xi-1)))
add("Li3(xi/(xi-1))", mp.polylog(3, xi/(xi-1)))
add("Li3(1/xi)", mp.polylog(3, 1/xi))
add("Li3((z-1)/(z-xi))", mp.polylog(3, (z-1)/(z-xi)))
add("Li2(z)L(1-z/xi)", mp.polylog(2, z)*L(1-z/xi))
add("Li2(z)L(1-z)", mp.polylog(2, z)*L(1-z))
add("Li2(z)L(z)", mp.polylog(2, z)*L(z))
add("Li2(z/xi)L(1-z/xi)", mp.polylog(2, z/xi)*L(1-z/xi))
add("Li2(z/xi)L(1-z)", mp.polylog(2, z/xi)*L(1-z))
add("Li2(z/xi)L(z)", mp.polylog(2, z/xi)*L(z))
add("Li2(z/xi)L(xi)", mp.polylog(2, z/xi)*L(xi))
add("L(z)L(1-z)L(1-z/xi)", L(z)*L(1-z)*L(1-z/xi))
add("L(z)^2 L(1-z/xi)", L(z)**2*L(1-z/xi))
add("L(z)L(1-z/xi)^2", L(z)*L(1-z/xi)**2)
add("L(1-z)^2 L(1-z/xi)", L(1-z)**2*L(1-z/xi))
add("L(1-z)L(1-z/xi)^2", L(1-z)*L(1-z/xi)**2)
add("L(1-z/xi)^3", L(1-z/xi)**3)
add("L(1-z/xi)^2 L(xi)", L(1-z/xi)**2*L(xi))
add("zeta3", mp.zeta(3))
add("zeta2 L(xi)", mp.pi**2/6*L(xi))
add("zeta2 L(1-z/xi)", mp.pi**2/6*L(1-z/xi))
add("zeta2 L(1-z)", mp.pi**2/6*L(1-z))
add("zeta2 L(z)", mp.pi**2/6*L(z))

vals = [v for _, v in cand]
rels = find_relations([Hz] + vals, maxcoeff=10**5, resid_decades=25, max_report=2)
print("relations found:", len(rels))
for row, mc, res in rels:
    print("maxcoeff", mc, "resid", mp.nstr(res, 2))
    terms = [(cand[i-1][0], row[i]) for i in range(1, len(row)) if row[i] != 0]
    print("H = ", [(n, c) for n, c in terms])
