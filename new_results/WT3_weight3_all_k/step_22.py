
import os
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK)
step2 = r'''"""WT3 step 2: exact W3(k) for k=3..15 from wt3_kernels.json, in the paper's weight-3 basis.
Controls: k=3 vs thm:M3 (term by term), k=4 vs rem:k4; quadrature at 70 dps for every k."""
import sys, os, json, time
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
sys.path.insert(0, WORK)
sys.path.insert(0, os.path.join(BASE, 'paper/squircle_note/scripts'))
from fractions import Fraction as Fr
import sympy as sp
import mpmath as mp
from qsqrt2 import QS2, SQ2, A, B
import polylog_exact as PE
from polylog_exact import P, LOG, LAM, LL, PI2, Z3, LI3, LI2
from w2_integrals import to_qs2
import kummer
mp.mp.dps = 70           # set AFTER all imports (dilog_exact import sets 60)
print("mp.mp.dps =", mp.mp.dps, flush=True)

ONE, ZERO = QS2(1), QS2(0)
v = sp.Symbol('v')
S2 = sp.sqrt(2)

data = json.load(open(os.path.join(WORK, 'wt3_kernels.json')))
def qsi_sp(pair):
    return sp.Rational(pair[0]) + sp.Rational(pair[1]) * S2
def kernel_sympy(k, key):
    terms = data['kernels'][str(k)][key]['real']
    return sum(qsi_sp(c) / (v - qsi_sp(p)) for p, c in terms)

# ---- exact log forms lt, l1t (as in m3_tails_exact) ----
LT = ({B: ONE, -A: ONE, ZERO: -ONE}, LAM * (-2))
L1T = ({A: ONE, -B: ONE, ZERO: -ONE}, LAM * (-2))
SUB = [(B, ONE), (ONE, A)]

def residues(R):
    res = {}
    for term in sp.Add.make_args(sp.expand(R)):
        num, den = sp.fraction(sp.together(term))
        dp = sp.Poly(sp.expand(den), v)
        assert dp.degree() == 1 and not sp.expand(num).has(v), term
        lc = dp.LC()
        p = to_qs2(sp.radsimp(-dp.nth(0) / lc))
        c = to_qs2(sp.radsimp(num / lc))
        res[p] = res.get(p, QS2(0)) + c
    back = sum(c.sympy() / (v - p.sympy()) for p, c in res.items())
    assert sp.simplify(sp.radsimp(sp.together(back - R))) == 0
    return res

def std_prod(l1, l2, z, v0, v1):
    (c1, k1), (c2, k2) = l1, l2
    tot = P()
    for x, cx in c1.items():
        for y, cy in c2.items():
            tot = tot + PE.int_loglog(x, y, z, v0, v1) * (cx * cy)
    for x, cx in c1.items():
        tot = tot + PE.int_log(x, z, v0, v1) * cx * k2
    for y, cy in c2.items():
        tot = tot + PE.int_log(y, z, v0, v1) * cy * k1
    tot = tot + PE.int_one(z, v0, v1) * k1 * k2
    return tot

def lam_at(lam, e):
    co, k = lam
    val = k
    for x, c in co.items():
        val = val + LOG(e - x) * c
    return val

def reg_prod(l1, l2, e, v0, v1):
    assert not lam_at(l1, e).d and not lam_at(l2, e).d
    tot = P()
    for x, cx in l1[0].items():
        for y, cy in l2[0].items():
            tot = tot + PE.int_reg2(x, y, e, v0, v1) * (cx * cy)
    return tot

def weight3_piece(i, j, R):
    l_a = LT
    l_b = LT if i == 2 else L1T
    assert (i, j) in ((2, 0), (1, 1))
    tot = P()
    for z, r in residues(R).items():
        for v0, v1 in SUB:
            if not (v0 <= z <= v1):
                tot = tot + std_prod(l_a, l_b, z, v0, v1) * r
            else:
                assert z in (v0, v1), z
                tot = tot + reg_prod(l_a, l_b, z, v0, v1) * r
    return tot

# numerics
r2 = mp.sqrt(2); aN, bN = 1 + r2, r2 - 1
lt = lambda x: mp.log((x - bN) * (x + aN) / (4 * x))
l1t = lambda x: mp.log((aN - x) * (x + bN) / (4 * x))
def quad_k(k):
    f20 = sp.lambdify(v, kernel_sympy(k, '20'), modules='mpmath')
    f11 = sp.lambdify(v, kernel_sympy(k, '11'), modules='mpmath')
    return mp.quad(lambda x: f20(x) * lt(x) ** 2 + f11(x) * lt(x) * l1t(x), [bN, 1, aN])

# ---- hardcoded thm:M3 W3' and rem:k4 ----
W3p_thm = (Z3 * QS2(Fr(-3, 16), Fr(7, 32)) + LI3(B) * QS2(1, 1) + LI3(-B) * QS2(1, -1)
           + LI3(QS2(0, Fr(1, 2))) * QS2(0, -2) + LI2(B) * (LL * 2 - LAM * QS2(0, 1))
           + PI2 * LL * QS2(Fr(-1, 8), Fr(1, 12)) + PI2 * LAM * QS2(Fr(-1, 24), Fr(-1, 24))
           + LAM * LAM * LAM * QS2(Fr(1, 24), Fr(1, 24)) + LAM * LAM * LL * QS2(0, Fr(1, 4))
           + LAM * LL * LL * QS2(Fr(-1, 2), Fr(-1, 2)) + LL * LL * LL * QS2(Fr(5, 6), Fr(-1, 6)))
DELTA_rem = Z3 * Fr(-1, 8) - PI2 * LAM * Fr(1, 24) - LAM * LAM * LAM * Fr(1, 12)

# Kummer relations once
KREL = [(f"Kummer (K) at (x,y)=({x}, {y})", kummer.kummer_relation(x, y)) for (x, y) in kummer.search()]
print(f"{len(KREL)} Kummer relations", flush=True)
basis = [B, -B, QS2(0, Fr(1, 2))]

def p_to_json(W):
    out = []
    for kkey, c in W.d.items():
        atoms = []
        for a in kkey:
            if isinstance(a, str):
                atoms.append(a)
            else:
                atoms.append([a[0], [str(a[1].p), str(a[1].q)]])
        out.append([atoms, [str(c.p), str(c.q)]])
    return out

values = {}
W3_of = {}
for k in range(3, 16):
    t0 = time.time()
    W3p = weight3_piece(2, 0, kernel_sympy(k, '20')) + weight3_piece(1, 1, kernel_sympy(k, '11'))
    q3 = quad_k(k)
    raws = W3p.num() - q3
    W3pb = PE.reduce_li2_all(W3p)
    W3pc, used, rk = PE.solve_li3_linear(W3pb, basis, extra_relations=KREL)
    resc = W3pc.num() - q3
    W3_of[k] = W3pc
    values[k] = p_to_json(W3pc)
    print(f"k={k}: |raw-quad|={mp.nstr(abs(raws),3)}  |clean-quad|={mp.nstr(abs(resc),3)}  "
          f"(dps={mp.mp.dps})  Li3-rel rank {rk}, {len(used)} used, {time.time()-t0:.1f}s", flush=True)

# controls
def pzero(W):
    return not W.d
c3 = pzero(W3_of[3] - W3p_thm)
print("CONTROL k=3 == thm:M3 W3' term by term:", "PASS" if c3 else "FAIL", flush=True)
c4 = pzero(W3_of[4] - (W3p_thm + DELTA_rem))
print("CONTROL k=4 == W3' - zeta3/8 - pi2 lam/24 - lam^3/12 (rem:k4):", "PASS" if c4 else "FAIL", flush=True)
DELTA = W3_of[4] - W3_of[3]
print("DELTA = W3(4)-W3(3) =", DELTA, flush=True)
print("DELTA == -zeta3/8 - pi2 lam/24 - lam^3/12 :", "PASS" if pzero(DELTA - DELTA_rem) else "FAIL", flush=True)
linok = all(pzero(W3_of[k] - (W3_of[3] + (W3_of[4] - W3_of[3]) * Fr(k - 3))) for k in range(3, 16))
print("LINEARITY W3(k) = W3(3)+(k-3)DELTA for k=3..15:", "PASS" if linok else "FAIL", flush=True)

with open(os.path.join(WORK, 'w3_values.json'), 'w') as fh:
    json.dump({"basis": ["z3", "Li3(b)", "Li3(-b)", "Li3(1/sqrt2)", "Li2(b)",
                          "lam=log2", "L=log(1+sqrt2)", "pi2"],
               "values": {str(k): values[k] for k in values}}, fh, indent=1)
print("saved w3_values.json", flush=True)
print("STEP2:", "PASS" if (c3 and c4 and linok) else "FAIL")
'''
with open('step_02.py','w') as fh:
    fh.write(step2)
print("written")
