"""step_30.py -- ladder derivation with real parts; then differentiation verification."""
import mpmath as mp
from ir_utils import find_relations
mp.mp.dps = 100

def Hval(z, xi):
    return mp.quad(lambda w: mp.polylog(2, w)/(w - xi), [0, z/2, z])

xi = mp.mpf('2.7'); z = mp.mpf('0.31')
Hz = Hval(z, xi)
L = mp.log
cand = [
 ("Li3(z)", mp.polylog(3, z)), ("Li3(z/xi)", mp.polylog(3, z/xi)),
 ("Li3(1-z)", mp.polylog(3, 1-z)), ("Li3((xi-z)/xi)", mp.polylog(3, (xi-z)/xi)),
 ("Li3(z/(z-xi))", mp.polylog(3, z/(z-xi))), ("Li3(xi/(xi-z))", mp.polylog(3, xi/(xi-z))),
 ("Li3((z-xi)/(z-1))", mp.polylog(3, (z-xi)/(z-1))), ("Li3((xi-z)/(xi-1))", mp.polylog(3, (xi-z)/(xi-1))),
 ("Li3(xi/(xi-1))", mp.polylog(3, xi/(xi-1))), ("Li3(1/xi)", mp.polylog(3, 1/xi)),
 ("Li3((z-1)/(z-xi))", mp.polylog(3, (z-1)/(z-xi))),
 ("Li2(z)L(1-z/xi)", mp.polylog(2, z)*L(1-z/xi)), ("Li2(z)L(1-z)", mp.polylog(2, z)*L(1-z)),
 ("Li2(z)L(z)", mp.polylog(2, z)*L(z)), ("Li2(z/xi)L(1-z/xi)", mp.polylog(2, z/xi)*L(1-z/xi)),
 ("Li2(z/xi)L(1-z)", mp.polylog(2, z/xi)*L(1-z)), ("Li2(z/xi)L(z)", mp.polylog(2, z/xi)*L(z)),
 ("Li2(z/xi)L(xi)", mp.polylog(2, z/xi)*L(xi)),
 ("LzL1zL1zxi", L(z)*L(1-z)*L(1-z/xi)), ("Lz2L1zxi", L(z)**2*L(1-z/xi)),
 ("LzL1zxi2", L(z)*L(1-z/xi)**2), ("L1z2L1zxi", L(1-z)**2*L(1-z/xi)),
 ("L1zL1zxi2", L(1-z)*L(1-z/xi)**2), ("L1zxi3", L(1-z/xi)**3),
 ("L1zxi2Lxi", L(1-z/xi)**2*L(xi)), ("zeta3", mp.zeta(3)),
 ("z2Lxi", mp.pi**2/6*L(xi)), ("z2L1zxi", mp.pi**2/6*L(1-z/xi)),
 ("z2L1z", mp.pi**2/6*L(1-z)), ("z2Lz", mp.pi**2/6*L(z)),
]
vals = [mp.re(v) for _, v in cand]
rels = find_relations([Hz] + vals, maxcoeff=10**5, resid_decades=20, max_report=2)
print("relations:", len(rels))
for row, mc, res in rels:
    print("maxcoeff", mc, "resid", mp.nstr(res, 2))
    print("H =", [(cand[i-1][0], row[i]) for i in range(1, len(row)) if row[i] != 0])
