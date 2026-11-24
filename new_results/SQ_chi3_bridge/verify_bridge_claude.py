"""Exact proof certificate for the chi3 bridge identity (Claude, 2026-10-22).

    48 chi3(b) = 32 sqrt2 L(chi8,3) - 3 pi^2 L + 4 L^3,   b = sqrt2 - 1, L = log(1 + sqrt2).

Both sides are written as Q-vectors on the 576 level-8 weight-3 H-words with the exact Cayley engine cmzv (A3f4/Q):
Li3(+-b) by cmzv.atom, L by atom(('log', A)), pi^2 by atom(('pi2',)), and
32 sqrt2 L(chi8,3) = 16[Li3(xi) + Li3(xi^7) - Li3(xi^3) - Li3(xi^5)] with depth-1 words H(-1,-1,j).
DIFF = LHS - RHS. Steps:
 1. numerics: each side against mpmath, and DIFF . values, with the 150-digit word values;
 2. membership mod two primes in the proved base (A3f3 weight3_with_lifts: double shuffle + distribution + weight-one +
    Rogers five-term lifts; and if needed the 13 Kummer rows over Q(sqrt2));
 3. EXACT certificate: DIFF = sum_i c_i row_i with rational c_i over an independent row subset, checked with Fractions on
    all 576 words.
Run from research/08_theory:  ../../.venv/bin/python SQ_chi3_bridge/verify_bridge_claude.py
"""
import sys, pickle, time
from fractions import Fraction as Fr
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
import mpmath as mp
import flint
from sympy.ntheory import sqrt_mod

mp.mp.dps = 150
X, Y = C.atom(('li', 3, C.B)), C.atom(('li', 3, -C.B))
LHS = C.add((X, 24), (Y, -24))
Hs = C.add(({(-1, -1, 1): Fr(1)}, 16), ({(-1, -1, 7): Fr(1)}, 16), ({(-1, -1, 3): Fr(1)}, -16), ({(-1, -1, 5): Fr(1)}, -16))
P2, Lv = C.atom(('pi2',)), C.atom(('log', C.A))
RHS = C.add((Hs, 1), (C.mul(P2, Lv), -3), (C.mul(Lv, C.mul(Lv, Lv)), 4))
DIFF = {w: c for w, c in C.add((LHS, 1), (RHS, -1)).items() if c != 0}

W = pickle.load(open('A3f3_level8_compute/weight3_with_lifts.pickle', 'rb'))
K = pickle.load(open('A3f3_level8_compute/k3_newHrows.pickle', 'rb'))['krows']
cols = list(W['cols']); idx = {c: i for i, c in enumerate(cols)}; n = len(cols)
assert all(w in idx for w in DIFF), "DIFF has a word outside the 576 columns"
print(f"DIFF: {len(DIFF)} words, max |coeff| {max(abs(c) for c in DIFF.values())}")

# 1. numerics
V = pickle.load(open('A3f6_level16/vals8_150digit.pickle', 'rb'))
val = lambda vec: mp.fsum(mp.mpf(c.numerator) / c.denominator * mp.mpc(*V[w]) for w, c in vec.items())
r2 = mp.sqrt(2); b = r2 - 1; L = mp.log(1 + r2)
chi3 = (mp.polylog(3, b) - mp.polylog(3, -b)) / 2
L8 = mp.nsum(lambda k: 1 / (8 * k + 1) ** 3 - 1 / (8 * k + 3) ** 3 - 1 / (8 * k + 5) ** 3 + 1 / (8 * k + 7) ** 3, [0, mp.inf])
print("LHS vec - 48 chi3(b)              :", mp.nstr(abs(val(LHS) - 48 * chi3), 3))
print("RHS vec - (32r2 L8 - 3pi^2 L + 4L^3):", mp.nstr(abs(val(RHS) - (32 * r2 * L8 - 3 * mp.pi ** 2 * L + 4 * L ** 3)), 3))
print("DIFF . values                     :", mp.nstr(abs(val(DIFF)), 3))


# 2. membership mod p (rational rows alone, then with Kummer rows)
def modp(x, p):
    x = Fr(x); return x.numerator % p * pow(x.denominator % p, -1, p) % p


def rows_modp(p, with_k):
    s = int(sqrt_mod(2, p))
    R = []
    for r in W['rows']:
        v = [0] * n
        for w, c in r: v[idx[w]] = (v[idx[w]] + modp(c, p)) % p
        R.append(v)
    if with_k:
        for r in K:
            v = [0] * n
            for w, (x, y) in r.items(): v[idx[w]] = (v[idx[w]] + modp(x, p) + s * modp(y, p)) % p
            R.append(v)
    return R


d = lambda p: [modp(DIFF.get(c, 0), p) for c in cols]
need_k = False
for p in (1000000033, 1000000241):
    for with_k in (False, True):
        R = rows_modp(p, with_k)
        r0 = flint.nmod_mat(len(R), n, [x for r in R for x in r], p).rank()
        r1 = flint.nmod_mat(len(R) + 1, n, [x for r in R for x in r] + d(p), p).rank()
        print(f"p={p} {'base+Kummer' if with_k else 'rational base'}: rank {r0}, +DIFF {r1}, in span: {r0 == r1}")
        if not with_k and r0 != r1: need_k = True

# 3. exact certificate over Q from the rational rows (if DIFF is in their span)
if need_k:
    print("DIFF needs the Kummer rows: exact Q(sqrt2) certificate not implemented here"); sys.exit(1)
p = 1000000033
R = rows_modp(p, False)
# independent rows: pivots of R^T
RT = flint.nmod_mat(n, len(R), [R[i][j] for j in range(n) for i in range(len(R))], p)
E, rk = RT.rref()
piv_rows = []
for i in range(rk):
    piv_rows.append(next(j for j in range(len(R)) if int(E[i, j]) != 0))
# pivot columns of the chosen rows
Msel = flint.nmod_mat(rk, n, [x for i in piv_rows for x in R[i]], p)
E2, rk2 = Msel.rref(); assert rk2 == rk
piv_cols = [next(j for j in range(n) if int(E2[i, j]) != 0) for i in range(rk)]
print(f"independent rows {rk}; solving {rk}x{rk} exactly over Q ...", flush=True)
t0 = time.time()
def accum(r):
    out = {}
    for w, a in r: out[w] = out.get(w, 0) + Fr(a)
    return out


rowsQ = [accum(W['rows'][i]) for i in piv_rows]
A = flint.fmpq_mat(rk, rk, [flint.fmpq(Fr(rowsQ[i].get(cols[j], 0)).numerator, Fr(rowsQ[i].get(cols[j], 0)).denominator)
                            for j in piv_cols for i in range(rk)])
bvec = flint.fmpq_mat(rk, 1, [flint.fmpq(Fr(DIFF.get(cols[j], 0)).numerator, Fr(DIFF.get(cols[j], 0)).denominator) for j in piv_cols])
sol = A.solve(bvec)
c = [Fr(int(sol[i, 0].p), int(sol[i, 0].q)) for i in range(rk)]
print(f"solved in {time.time() - t0:.0f}s; {sum(1 for x in c if x)} nonzero multipliers")
# exact check on ALL 576 words
comb = {}
for ci, r in zip(c, rowsQ):
    if ci:
        for w, a in r.items(): comb[w] = comb.get(w, 0) + ci * Fr(a)
comb = {w: a for w, a in comb.items() if a != 0}
ok = comb == DIFF
print("EXACT CERTIFICATE:", "DIFF = sum c_i row_i holds on all 576 words" if ok else "FAILED")
pickle.dump({'piv_rows': piv_rows, 'multipliers': c, 'DIFF': DIFF, 'ok': ok},
            open('SQ_chi3_bridge/bridge_certificate.pickle', 'wb'))
