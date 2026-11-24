
import json
from mpmath import mp, mpf, zeta, fabs, log, matrix, lu_solve, nstr, pslq
mp.dps = 125
exec(open('av2core.py').read())
hifi = json.load(open('b3_hifi.json'))
b2tab = json.load(open('b2_table.json'))
z3, z2 = zeta(3), zeta(2)
ps = [24, 48, 64, 256, 1024, 4096, 16384]
for p in [16384]:
    b2tab[str(p)] = dict(B2 = mp.nstr(B2_quad2(p), 115))
E = {}
print("=== E(p) = B3(p) - 2 B2(p) ===")
for p in ps:
    B3 = mpf(hifi[str(p)]['B3'])
    B2 = mpf(b2tab[str(p)]['B2'])
    E[p] = B3 - 2*B2
    print("p=%6d  E=%s   E*p^2=%s" % (p, nstr(E[p], 16), nstr(E[p]*p**2, 14)))
# model: E = (e2 + f2 ln p) e^2 + (e3 + f3 ln p) e^3 + (e4 + f4 ln p) e^4,   e=1/p
# staged: fit on 4 largest with 6 unknowns, validate on rest
def fit(pvals, nunk):
    A = []; b = []
    for p in pvals:
        L = log(p); e = mpf(1)/p
        row = [e**2, L*e**2, e**3, L*e**3]
        if nunk >= 6: row += [e**4, L*e**4]
        A.append(row); b.append(E[p])
    return lu_solve(matrix(A), matrix(b))
sol = fit([1024, 4096, 16384, 256], 4)
e2, f2, e3, f3 = sol
print("staged fit (256..16384):")
print("e2 =", nstr(e2, 40))
print("f2 =", nstr(f2, 40), "  vs -2 zeta(2):", nstr(f2+2*z2, 10))
print("e3 =", nstr(e3, 40), "  vs -5 zeta(3):", nstr(e3+5*z3, 10), " vs -zeta(3):", nstr(e3+z3,10))
print("f3 =", nstr(f3, 40))
def model(p, s):
    L = log(p); e = mpf(1)/p
    return (s[0]+s[1]*L)*e**2 + (s[2]+s[3]*L)*e**3
for p in [24, 48, 64]:
    print("validate p=%4d: model-E = %s" % (p, nstr(model(p, sol)-E[p], 8)))
# A2 = zeta(2) - e2  (from B3 = 2 zeta3 + [2 zeta2 ln p - A2] e^2 ...)
print("=> A2 = zeta2 - e2 =", nstr(z2 - e2, 40))
print("=== PSLQ on e2 (and candidates) ===")
basis = [mpf(1), z2, z3, z2**2, z2*z3, z3**2, zeta(4), zeta(5), log(2), z2*log(2), z3*log(2), log(2)**2]
for nm, val in [('e2', e2), ('e3', e3), ('A2', z2-e2)]:
    rel = pslq([val]+basis, tol=mpf('1e-28'), maxcoeff=10**10, maxsteps=2000)
    print(nm, "=", nstr(val, 32), "->", rel)
