
import os, json, sys
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK)
sys.path.insert(0, WORK); sys.path.insert(0, os.path.join(BASE,'paper/squircle_note/scripts'))
from fractions import Fraction as Fr
import sympy as sp
import mpmath as mp
import fastred
from buildintegrand import build_integrand_pf
from qsqrt2 import QS2, SQ2, A, B
import polylog_exact as PE
from polylog_exact import P, LOG, LAM, LL, PI2, Z3, LI3, LI2
from w2_integrals import to_qs2
import kummer
mp.mp.dps = 70
ONE, ZERO = QS2(1), QS2(0)
v = sp.Symbol('v'); S2 = sp.sqrt(2)

red = fastred.reduce_k_fast(2, integrand_pf=build_integrand_pf(2))
R20 = fastred.pf_to_sympy(red[(2,0)]); R11 = fastred.pf_to_sympy(red[(1,1)])

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
        p = to_qs2(sp.radsimp(-dp.nth(0) / lc)); c = to_qs2(sp.radsimp(num / lc))
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
    l_a = LT; l_b = LT if i == 2 else L1T
    tot = P()
    for z, r in residues(R).items():
        for v0, v1 in SUB:
            if not (v0 <= z <= v1):
                tot = tot + std_prod(l_a, l_b, z, v0, v1) * r
            else:
                assert z in (v0, v1)
                tot = tot + reg_prod(l_a, l_b, z, v0, v1) * r
    return tot

r2 = mp.sqrt(2); aN, bN = 1 + r2, r2 - 1
lt = lambda x: mp.log((x - bN) * (x + aN) / (4 * x))
l1t = lambda x: mp.log((aN - x) * (x + bN) / (4 * x))
f20 = sp.lambdify(v, R20, modules='mpmath'); f11 = sp.lambdify(v, R11, modules='mpmath')
q = mp.quad(lambda x: f20(x) * lt(x)**2 + f11(x) * lt(x) * l1t(x), [bN, 1, aN])

W3p = weight3_piece(2, 0, R20) + weight3_piece(1, 1, R11)
print("|raw-quad| =", mp.nstr(abs(W3p.num() - q), 3), "dps", mp.mp.dps)
KREL = [(f"K({x},{y})", kummer.kummer_relation(x, y)) for (x, y) in kummer.search()]
W3pc, used, rk = PE.solve_li3_linear(PE.reduce_li2_all(W3p), [B, -B, QS2(0, Fr(1,2))], extra_relations=KREL)
print("W3(2) =", W3pc)
print("|clean-quad| =", mp.nstr(abs(W3pc.num() - q), 3))

W3_M2 = (Z3 * QS2(Fr(-19,16), Fr(189,64)) + LI3(B) * QS2(11, 15) + LI3(-B) * QS2(11, -15)
         + LI3(QS2(0, Fr(1,2))) * QS2(0, -27) + LI2(B) * (LL * 12 - LAM * QS2(0, Fr(7,2)))
         + PI2 * LL * QS2(Fr(-7,12), Fr(15,8)) + PI2 * LAM * QS2(Fr(-1,24), Fr(-65,48))
         + LAM**3 * QS2(Fr(1,6), Fr(9,16)) + LAM**2 * LL * QS2(0, Fr(11,8))
         + LAM * LL**2 * QS2(-1, Fr(-7,4)) + LL**3 * QS2(Fr(25,6), Fr(-7,3)))
diff = W3pc - W3_M2
print("W3(2) == thm:M2 weight-3 part:", "PASS" if not diff.d else ("FAIL: " + str(diff)))

vals = json.load(open('w3_values.json'))['values']
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
W3_3 = p_from_json(vals['3'])
print("W3(3) - W3(2) =", W3_3 - W3pc)
print("(for k>=3 the step is constant = -1/8*z3 - 1/24*lam*pi2 - 1/12*lam^3)")
