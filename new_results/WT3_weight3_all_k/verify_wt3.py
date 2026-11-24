"""WT3 verification: reruns tasks 1-3 (+ the task-4 reduction identities) and prints PASS/FAIL.
Target < 15 min. Prints mp.mp.dps next to every numerical residual."""
import sys, os, json, time
BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
sys.path.insert(0, WORK)
sys.path.insert(0, os.path.join(BASE, 'paper/squircle_note/scripts'))
from fractions import Fraction as Fr
import sympy as sp
import mpmath as mp
import fastred
from fastred import (QSI, ZERO, ONE, IQ, AQ, BQ, pf_zero, pf_add, pf_sub, pf_scale, pf_split,
                     pf_mul_DL, pf_eval_at, pf_from_sympy, pf_to_sympy, pf_fold, DL1, DL2, pf_neg)  # pf_neg added (M1 session)
from buildintegrand import build_integrand_pf
from qsqrt2 import QS2, A, B
import polylog_exact as PE
from polylog_exact import P, LOG, LAM, LL, PI2, Z3, LI3, LI2
from w2_integrals import to_qs2
import kummer
mp.mp.dps = 65    # set AFTER imports (dilog_exact import sets 60)
T0 = time.time()
results = []
def report(name, ok, extra=""):
    results.append((name, ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  {extra}  ({time.time()-T0:.0f}s)", flush=True)

def pf_is_zero(P):
    return not P[0] and all(e.is_zero() for e in P[1]) and all(e.is_zero() for e in P[2]) and P[3].is_zero()

# ---------- task 1: kernels ----------
kernels = {}
for k in range(3, 16):
    red = fastred.reduce_k_fast(k, integrand_pf=build_integrand_pf(k))
    kernels[k] = {key: red.get(key, pf_zero()) for key in [(2,0),(1,1)]}
    for key, PFk in kernels[k].items():  # renamed (M1 session, 2026-10-27): 'P' shadowed polylog_exact.P
        assert not PFk[1] and not PFk[2] and PFk[3].is_zero()
# saved json agrees
data = json.load(open(os.path.join(WORK,'wt3_kernels.json')))
ok = True
for k in range(3, 16):
    for key, nm in [((2,0),'20'), ((1,1),'11')]:
        saved = pf_zero()
        for p, c in data['kernels'][str(k)][nm]['real']:
            rq = QSI(Fr(p[0]), Fr(p[1])); cq = QSI(Fr(c[0]), Fr(c[1]))
            saved[0][rq] = cq
        ok &= pf_is_zero(pf_sub(saved, kernels[k][key]))
report("task1: kernels k=3..15 == wt3_kernels.json", ok)
# affine pattern
S20 = pf_zero(); S20[0][-AQ] = QSI(Fr(1,2)); S20[0][-BQ] = QSI(Fr(-1,2))
S11 = pf_zero(); S11[0][ZERO] = QSI(Fr(1,2)); S11[0][-AQ] = QSI(-1)
ok = all(pf_is_zero(pf_sub(kernels[k][(2,0)], pf_add(kernels[3][(2,0)], pf_scale(S20, Fr(k-3)))))
         and pf_is_zero(pf_sub(kernels[k][(1,1)], pf_add(kernels[3][(1,1)], pf_scale(S11, Fr(k-3)))))
         for k in range(3, 16))
report("task3(b): kernels affine in k for k>=3 (exact)", ok)
# control vs slow ibp_reduce at k=3,4 (degree-2 keys)
from ibp_reduce import reduce_k
ok = True
for k in [3, 4]:
    slow = reduce_k(k, verbose=False)
    for key in [(2,0),(1,1)]:
        ok &= pf_is_zero(pf_sub(fastred.pf_from_sympy(slow.get(key, 0)), kernels[k][key]))
report("task1 control: slow ibp_reduce == fast reducer, degree-2 keys k=3,4", ok)

# ---------- task 2: values ----------
ONEQ, ZERQ = QS2(1), QS2(0)
LT = ({B: ONEQ, -A: ONEQ, ZERQ: -ONEQ}, LAM * (-2))
L1T = ({A: ONEQ, -B: ONEQ, ZERQ: -ONEQ}, LAM * (-2))
SUB = [(B, ONEQ), (ONEQ, A)]
v = fastred.v
def residues(R):
    res = {}
    for term in sp.Add.make_args(sp.expand(R)):
        num, den = sp.fraction(sp.together(term))
        dp = sp.Poly(sp.expand(den), v)
        lc = dp.LC()
        p = to_qs2(sp.radsimp(-dp.nth(0)/lc)); c = to_qs2(sp.radsimp(num/lc))
        res[p] = res.get(p, QS2(0)) + c
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
    return tot + PE.int_one(z, v0, v1) * k1 * k2
def lam_at(lam, e):
    co, kk = lam
    val = kk
    for x, c in co.items():
        val = val + LOG(e - x) * c
    return val
def reg_prod(l1, l2, e, v0, v1):
    tot = P()
    for x, cx in l1[0].items():
        for y, cy in l2[0].items():
            tot = tot + PE.int_reg2(x, y, e, v0, v1) * (cx * cy)
    return tot
def w3val(k):
    tot = P()
    for (i, j), l_b in [((2,0), LT), ((1,1), L1T)]:
        R = pf_to_sympy(kernels[k][(i,j)])
        for z, r in residues(R).items():
            for v0, v1 in SUB:
                if not (v0 <= z <= v1):
                    tot = tot + std_prod(LT, l_b, z, v0, v1) * r
                else:
                    tot = tot + reg_prod(LT, l_b, z, v0, v1) * r
    return tot
r2 = mp.sqrt(2); aN, bN = 1+r2, r2-1
lt = lambda x: mp.log((x-bN)*(x+aN)/(4*x)); l1t = lambda x: mp.log((aN-x)*(x+bN)/(4*x))
def quad_k(k):
    f20 = sp.lambdify(v, pf_to_sympy(kernels[k][(2,0)]), modules='mpmath')
    f11 = sp.lambdify(v, pf_to_sympy(kernels[k][(1,1)]), modules='mpmath')
    return mp.quad(lambda x: f20(x)*lt(x)**2 + f11(x)*lt(x)*l1t(x), [bN, 1, aN])
KREL = [(f"K({x},{y})", kummer.kummer_relation(x, y)) for (x, y) in kummer.search()]
basis = [B, -B, QS2(0, Fr(1,2))]
W3p_thm = (Z3 * QS2(Fr(-3,16), Fr(7,32)) + LI3(B)*QS2(1,1) + LI3(-B)*QS2(1,-1)
           + LI3(QS2(0,Fr(1,2)))*QS2(0,-2) + LI2(B)*(LL*2 - LAM*QS2(0,1))
           + PI2*LL*QS2(Fr(-1,8), Fr(1,12)) + PI2*LAM*QS2(Fr(-1,24),Fr(-1,24))
           + LAM*LAM*LAM*QS2(Fr(1,24),Fr(1,24)) + LAM*LAM*LL*QS2(0,Fr(1,4))
           + LAM*LL*LL*QS2(Fr(-1,2),Fr(-1,2)) + LL*LL*LL*QS2(Fr(5,6),Fr(-1,6)))
DELTA_rem = Z3*Fr(-1,8) - PI2*LAM*Fr(1,24) - LAM*LAM*LAM*Fr(1,12)
W3 = {}
quadres = []
for k in range(3, 16):
    W3p = w3val(k)
    q = quad_k(k)
    W3pc, used, rk = PE.solve_li3_linear(PE.reduce_li2_all(W3p), basis, extra_relations=KREL)
    W3[k] = W3pc
    quadres.append(abs(W3pc.num() - q))
maxres = max(quadres)
report("task2: W3(k) exact vs quadrature, all k=3..15", maxres < mp.mpf(10)**-60,
       f"max residual {mp.nstr(maxres,3)} at dps={mp.mp.dps}")
report("task2 control k=3 == thm:M3 W3'", not (W3[3]-W3p_thm).d)
report("task2 control k=4 == rem:k4", not (W3[4]-(W3p_thm+DELTA_rem)).d)
DELTA = W3[4]-W3[3]
ok = all(not (W3[k]-(W3[3]+DELTA*Fr(k-3))).d for k in range(3,16))
report("task3: W3(k) = W3(3)+(k-3)DELTA, DELTA=-z3/8-pi2 lam/24-lam^3/12 (exact)", ok and not (DELTA-DELTA_rem).d)
# saved json agrees
vals = json.load(open(os.path.join(WORK,'w3_values.json')))['values']
def p_from_json(lst):
    out = P()
    for atoms, c in lst:
        term = P.const(QS2(Fr(c[0]), Fr(c[1])))
        for a in atoms:
            if isinstance(a, str):
                term = term * P({(a,): 1})
            else:
                z = QS2(Fr(a[1][0]), Fr(a[1][1]))
                term = term * (LI2(z) if a[0]=='Li2' else LI3(z))
        out = out + term
    return out
ok = all(not (W3[k]-p_from_json(vals[str(k)])).d for k in range(3,16))
report("task2: values == w3_values.json", ok)

# ---------- task 4 reduction: W1, W2, C1-C4 ----------
record_ok = True
for k in range(3, 8):   # fixed-point degree-2 equations on a few levels (fast)
    red = {}
    pf_int = build_integrand_pf(k)
    work = {}
    for key, R in pf_int.items():
        kk, RR = fastred.fold_key(key, R)
        work[kk] = pf_add(work.get(kk, pf_zero()), RR)
    gms = {}
    Ws = {}
    for deg in range(k, -1, -1):
        for key in sorted([kk for kk in work if sum(kk)==deg], reverse=True):
            R = work.pop(key)
            if deg == 0: continue
            simple, higher, const = pf_split(R)
            Ws[key] = R
            if not higher[1] and not higher[2]: continue
            G = fastred.pf_integrate_higher(higher)
            C = pf_eval_at(G, BQ)
            Gm = [G[0], G[1], G[2], -C]
            gms[key] = Gm
            i, j = key
            if i >= 1:
                kk, RR = fastred.fold_key((i-1,j), pf_scale(pf_mul_DL(Gm, DL1), Fr(-i)))
                work[kk] = pf_add(work.get(kk, pf_zero()), RR)
            if j >= 1:
                kk, RR = fastred.fold_key((i,j-1), pf_scale(pf_mul_DL(Gm, DL2), Fr(-j)))
                work[kk] = pf_add(work.get(kk, pf_zero()), RR)
    Gm30, Gm21 = gms[(3,0)], gms[(2,1)]
    rhs20 = pf_add(pf_scale(pf_mul_DL(Gm30, DL1), Fr(-3)), pf_neg(pf_mul_DL(Gm21, DL2)))
    rhs11 = pf_scale(pf_mul_DL(Gm21, DL1), Fr(-2))
    record_ok &= pf_is_zero(pf_sub(Ws[(2,0)], rhs20)) and pf_is_zero(pf_sub(Ws[(1,1)], rhs11))
    if k == 3:
        G3, G21 = Gm30, Gm21
report("task4 (P2): degree-2 fixed-point equations k=3..7", record_ok)
s30 = [{}, [QSI(0,Fr(1,12))], [QSI(0,Fr(1,12))], QSI(Fr(-1,12))]
s21 = [{}, [QSI(0,Fr(-1,4))], [QSI(0,Fr(-1,4))], QSI(Fr(1,4))]
R20_3 = pf_zero(); R20_3[0][ZERO] = QSI(Fr(1,4),Fr(-1,4)); R20_3[0][AQ] = QSI(Fr(1,4)); R20_3[0][-AQ] = QSI(Fr(-1,4)); R20_3[0][-BQ] = QSI(Fr(-1,2))
R11_3 = pf_zero(); R11_3[0][ZERO] = QSI(Fr(1,4),Fr(1,4)); R11_3[0][-AQ] = QSI(Fr(-1,2))
def simple_part(P):
    s, h, c = pf_split(P)
    return s
C1 = pf_is_zero(pf_sub(simple_part(pf_add(pf_scale(pf_mul_DL(G3, DL1), Fr(-3)), pf_neg(pf_mul_DL(G21, DL2)))), R20_3))  # G21 = Gm21^(3) (M1 session fix: was the stale k=7 Gm21)
C2 = pf_is_zero(pf_sub(simple_part(pf_add(pf_scale(pf_mul_DL(s30, DL1), Fr(-3)), pf_neg(pf_mul_DL(s21, DL2)))), S20))
C3 = pf_is_zero(pf_sub(simple_part(pf_scale(pf_mul_DL(G21, DL1), Fr(-2))), R11_3))
C4 = pf_is_zero(pf_sub(simple_part(pf_scale(pf_mul_DL(s21, DL1), Fr(-2))), S11))
report("task4 (P3): reduction identities C1-C4 (exact)", C1 and C2 and C3 and C4)
print(f"\nverify_wt3: {sum(1 for _,ok in results if ok)}/{len(results)} checks PASS, {time.time()-T0:.0f}s total, dps={mp.mp.dps}")
bad = [n for n, ok in results if not ok]
print("OVERALL:", "PASS" if not bad else f"FAIL: {bad}")
