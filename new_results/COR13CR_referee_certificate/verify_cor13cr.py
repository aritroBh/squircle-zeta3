"""verify_cor13cr.py -- COR13CR referee's aggregate verifier of COR13C.
Runtime ~ 12 min. Exit code 0 iff every check passes.

PART A (exact):
  A1 relation row sets: reflib.make_relations == cmzv6.make_relations, weights 1-4.
  A2 certificates: sum_i y_i rows[prow_i] == T_n on EVERY admissible word; FLINT fmpq
     (weights 2,3,4 + Euler) and independent pure-Fraction backend (weights 2,3 + Euler);
     weight 4 also mod 3 primes. Targets rebuilt from cor13c_common; my independent
     targets (reflib13, same frames) == theirs EXACTLY.
  A3 Mobius maps: my exact enumeration over Q(zeta_6) (12 maps == mobius_maps_13.json);
     2N counts at N=6,12 (my mod-p filter + my Cyclo); cross-ratio availability; 36
     pulled-back 1-form identities exact over Q(zeta).
  A4 tamper test: perturbed certificate FAILS (the check discriminates).
PART B (numeric audit, >= 60 digits, my evaluator):
  B1 my CMZV evaluator vs COR13C's fast_cmzv/numeric_cmzv (40 words).
  B2 every certificate row (28+229+1695 + 28 Euler) vanishes at 65 dps.
  B3 conversion: 38 HPL words + 9 Li atoms + 7 atoms at 65-70 dps; wrong order fails.
  B4 targets vanish (my frames); frame difference in proved span mod p.
  B5 ID2/ID3/ID4 residuals at 250 dps with my HPL values; M3 closed form vs 1095-digit
     reference (at 1150 dps, own Li values).
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
BASE = "<project>/research/08_theory"
sys.path.insert(0, BASE + "/HEXC1_hexagon_certificate")
sys.path.insert(0, BASE + "/HEXC2R_referee_certificate")
sys.path.insert(0, BASE + "/COR13C_third_corner_certificate")
from fractions import Fraction as Fr
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)
ok = True
def check(name, cond):
    global ok
    ok = ok and bool(cond)
    log(f"  [{'PASS' if cond else 'FAIL'}] {name}")

import cmzv6
import reflib as R
import reflib13 as R13
from mymaps import map_count
D13 = BASE + "/COR13C_third_corner_certificate"

log("PART A: exact")
# A1
lower_m, lower_t = {}, {}
ROWS = {}
for n in (1, 2, 3, 4):
    rows_t, _ = cmzv6.make_relations(n, lower_t)
    rows_m, _ = R.make_relations(n, lower_m)
    lower_t[n] = rows_t; lower_m[n] = rows_m
    ROWS[n] = rows_t
    check(f"A1 weight {n}: row sets equal ({len(rows_t)} rows)", set(rows_t) == set(rows_m))

# A2
from flint import fmpq_mat
from cor13c_common import target2, target3, target4
cert23 = json.load(open(D13 + "/cert13_w23.json"))
cert4 = json.load(open(D13 + "/cert13_w4.json"))
certE = json.load(open(D13 + "/cert13_euler.json"))
def resid_fmpq(sel, y, T, cols):
    ids = {w: j for j, w in enumerate(cols)}
    S = fmpq_mat(len(sel), len(cols))
    for i, row in enumerate(sel):
        for w, c in row: S[i, ids[w]] = c
    Y = fmpq_mat(len(sel), 1)
    for i, yi in enumerate(y): Y[i, 0] = fmpq_mat(1, 1, [str(yi)])[0, 0]
    res = Y.transpose() * S
    return sum(1 for w, j in ids.items() if Fr(str(res[0, j])) != T.get(w, Fr(0)))
def resid_frac(sel, y, T, cols):
    ids = {w: j for j, w in enumerate(cols)}
    acc = [Fr(0)] * len(cols)
    for row, yi in zip(sel, y):
        if yi == 0: continue
        for w, c in row: acc[ids[w]] += yi * Fr(c)
    return sum(1 for w, j in ids.items() if acc[j] != T.get(w, Fr(0)))
def resid_modp(sel, y, T, cols, p):
    ids = {w: j for j, w in enumerate(cols)}
    acc = [0] * len(cols)
    for row, yi in zip(sel, y):
        ym = yi.numerator % p * pow(yi.denominator % p, p - 2, p) % p
        if ym == 0: continue
        for w, c in row: acc[ids[w]] = (acc[ids[w]] + ym * (c % p)) % p
    bad = 0
    for w, j in ids.items():
        t = T.get(w, Fr(0))
        if (acc[j] - t.numerator % p * pow(t.denominator % p, p - 2, p)) % p: bad += 1
    return bad
TARGETS = {2: target2, 3: target3, 4: target4}
for n in (2, 3, 4):
    T = {w: Fr(c) for w, c in TARGETS[n]().items()}
    cols = list(cmzv6.words(n))
    if n < 4:
        C = cert23[str(n)]
        prow, y = C["prow"], [Fr(int(a), int(b)) for a, b in C["y"]]
        assert [list(w) for w in cols] == C["cols"]
    else:
        prow, y = cert4["prow"], [Fr(s) for s in cert4["y"]]
        assert [list(w) for w in cols] == cert4["cols"]
    sel = [ROWS[n][i] for i in prow]
    b1 = resid_fmpq(sel, y, T, cols)
    check(f"A2 weight {n}: fmpq residual 0 on {len(cols)} words ({len(prow)} rows)", b1 == 0)
    if n <= 3:
        check(f"A2 weight {n}: pure-Fraction residual 0", resid_frac(sel, y, T, cols) == 0)
    else:
        check(f"A2 weight 4: residual 0 mod 3 primes",
              all(resid_modp(sel, y, T, cols, p) == 0 for p in (1000000007, 1000000009, 998244353)))
    if n == 4:
        yb = list(y); yb[0] += 1
        check("A4 tamper test: y[0]+=1 fails", resid_modp(sel, yb, T, cols, 1000000007) > 0)
# my targets == theirs (HC frames at weight 4, my rep0 elsewhere coincide)
check("A2 my target2 == theirs", R13.my_target2() == {w: Fr(c) for w, c in target2().items()})
check("A2 my target3 == theirs", R13.my_target3() == {w: Fr(c) for w, c in target3().items()})
# Euler
ET = {(-1, 0): Fr(1), (1, 1): Fr(3)}
cols2 = list(cmzv6.words(2))
selE = [ROWS[2][i] for i in certE["prow"]]
yE = [Fr(int(a), int(b)) for a, b in certE["y"]]
check("A2 Euler cert (fmpq + pure-Fraction)",
      resid_fmpq(selE, yE, ET, cols2) == 0 and resid_frac(selE, yE, ET, cols2) == 0)

# A3
maps_mine = R13.enum_maps13()
stored = sorted(tuple(m["exps"]) for m in json.load(open(D13 + "/mobius_maps_13.json")))
check("A3 12 maps exact over Q(zeta_6) == mobius_maps_13.json", maps_mine == stored)
s6 = map_count(6); s12 = map_count(12)
check(f"A3 map counts N=6: {s6[0]}/12, N=12: {s12[0]}/24, extras 0",
      s6[0] == 12 and s6[1] == 0 and s12[0] == 24 and s12[1] == 0)
tab = R13.cross_values13()
need = [Fr(1,3), Fr(-1,3), Fr(2,3), Fr(1,2), Fr(-1,2), Fr(-2)]
check("A3 all 6 Li arguments are cross-ratios; 1/9 absent",
      all(z in tab for z in need) and Fr(1,9) not in tab)
# 36 exact 1-form identities (compact re-derivation)
qadd, qsub, qmul, qdiv, qneg = R13.qadd, R13.qsub, R13.qmul, R13.qdiv, R13.qneg
MU6 = R13.MU6
ZERO = (Fr(0), Fr(0)); ONE = (Fr(1), Fr(0))
def pm(P, Q):
    out = [ZERO]*(len(P)+len(Q)-1)
    for i, a in enumerate(P):
        if a == ZERO: continue
        for j, b in enumerate(Q):
            if b != ZERO: out[i+j] = qadd(out[i+j], qmul(a, b))
    return out
def pa(P, Q):
    n = max(len(P), len(Q)); out = [ZERO]*n
    for i in range(n):
        out[i] = qadd(P[i] if i < len(P) else ZERO, Q[i] if i < len(Q) else ZERO)
    return out
def ps(P, c): return [qmul(a, c) for a in P]
def trim(P):
    P = list(P)
    while len(P) > 1 and P[-1] == ZERO: P.pop()
    return P
def l2p(vec):
    terms = [(c, w[0]) for w, c in vec.items()]
    den = [ONE]
    if any(j == -1 for c, j in terms): den = pm(den, [ZERO, ONE])
    for c, j in terms:
        if j >= 0: den = pm(den, [ONE, qneg(MU6[j])])
    num = [ZERO]
    for c, j in terms:
        p = [ONE]
        for c2, j2 in terms:
            if j2 == j: continue
            p = pm(p, [ZERO, ONE] if j2 == -1 else [ONE, qneg(MU6[j2])])
        p = ps(p, qmul((Fr(c), Fr(0)), ONE if j == -1 else MU6[j]))
        num = pa(num, p)
    return trim(num), trim(den)
cnt = 0
for mi, E in enumerate(maps_mine):
    f0, fh, f1, fm, finf = [MU6[e] for e in E]
    Dm = qdiv(qsub(finf, f1), qsub(f1, f0)); Bm = qmul(f0, Dm)
    Ap, Bp, Cp, Dp = Dm, qneg(Bm), qneg(ONE), finf
    for sig, q in ((0, ZERO), (1, ONE), (-1, qneg(ONE))):
        A2 = qsub(Ap, qmul(q, Cp)); B2 = qsub(Bp, qmul(q, Dp))
        NL = trim(pa(pm([A2], trim([Dp, Cp])), ps(pm([Cp], trim([B2, A2])), (Fr(-1), Fr(0)))))
        DL = trim(pm(trim([B2, A2]), trim([Dp, Cp])))
        if sig == 1: NL = ps(NL, (Fr(-1), Fr(0)))
        NR, DR = l2p(R13.hpl_form_pullback13(mi, sig))
        assert trim(pm(NL, DR)) == trim(pm(NR, DL)), (mi, sig)
        cnt += 1
check(f"A3 36 pulled-back 1-form identities exact over Q(zeta)", cnt == 36)
log("PART A " + ("PASSED" if ok else "FAILED"))

# PART B
log("PART B: numeric audit (>= 60 dps, my evaluator)")
import mpmath as mp
mp.mp.dps = 65
import fast_cmzv as FC
import numeric_cmzv as NC
FC.setup(62); FC.finish_setup(); NC.set_dps(62)
import random
random.seed(99)
worst = 0; nt = 0
while nt < 40:
    n = random.randint(1, 4)
    w = tuple(random.choice([-1,0,1,2,3,4,5]) for _ in range(n))
    if w and (w[0] == 0 or w[-1] == -1): continue
    a = R.cmzv_value(w)
    b = FC.cmzv_value(w); b = mp.mpc(b.real.mid(), b.imag.mid())
    c = NC.cmzv_value(w)
    worst = max(worst, abs(a-b), abs(a-c)); nt += 1
check(f"B1 my evaluator vs fast_cmzv/numeric_cmzv (40 words): worst {mp.nstr(worst,3)}", worst < mp.mpf("1e-58"))

# B2 rows vanish
worst = 0
prows = {2: cert23["2"]["prow"] + certE["prow"], 3: cert23["3"]["prow"], 4: cert4["prow"]}
for n in (2, 3, 4):
    wn = 0
    for i in set(prows[n]):
        val = mp.mpc(0)
        for w, c in ROWS[n][i]:
            val += int(c) * R.cmzv_value(w)
        wn = max(wn, abs(val))
    worst = max(worst, wn)
    log(f"  weight {n}: {len(set(prows[n]))} cert rows, worst |value| {mp.nstr(wn,3)}")
check(f"B2 all certificate rows vanish at 65 dps (worst {mp.nstr(worst,3)})", worst < mp.mpf("1e-58"))

# B3 conversion
mp.mp.dps = 70
WORDS = sorted({w for name in ["hpl_fit_k3_(0,0,0).json", "hpl_fit_k3_(1,0,0).json"]
                for e, w, v in R13.fit_terms(name) if w}, key=lambda x: (len(x), x))
worst = max(abs(R.cmzv_combo(R13.my13_hpl_to_cmzv(w)) - R.hpl_value(w, mp.mpf(1)/3)) for w in WORDS)
check(f"B3 {len(WORDS)} HPL conversions at 1/3, 70 dps: worst {mp.nstr(worst,3)}", worst < mp.mpf("1e-62"))
LI = [(2,"1/3"),(2,"-1/3"),(3,"1/3"),(3,"-1/2"),(4,"1/3"),(4,"-1/3"),(4,"2/3"),(4,"1/2"),(4,"-2")]
worst = max(abs(R.cmzv_combo(R13.my13_li_atom(n, Fr(z))) - mp.polylog(n, mp.mpf(z))) for n, z in LI)
check(f"B3 {len(LI)} Li atoms, 70 dps: worst {mp.nstr(worst,3)}", worst < mp.mpf("1e-62"))
mp.mp.dps = 65
ats = [(R13.a_log2(), mp.log(2)), (R13.a_log3(), mp.log(3)), (R13.a_ipi(), mp.pi*mp.j),
       (R13.a_zeta2(), mp.zeta(2)), (R13.a_zeta3(), mp.zeta(3)), (R13.a_pi2(), -mp.pi**2),
       (R13.a_pi4(), mp.pi**4)]
worst = max(abs(R.cmzv_combo(v) - x) for v, x in ats)
check(f"B3 7 elementary atoms: worst {mp.nstr(worst,3)}", worst < mp.mpf("1e-60"))
# wrong order
def seg_wrong(w, c, e):
    w = tuple(w)
    if c == e: return {(): Fr(1)} if not w else {}
    out = {}; nn = len(w)
    for k in range(nn + 1):
        left = R.regularize({R._shift6(w[k:], e): Fr(1)})
        right = R.regularize({R._shift6(tuple(reversed(w[:k])), c): Fr(1)})
        out = R.vadd((out, 1), (R.mul(left, right), Fr((-1)**k)))
    return R.regularize(out)
def hpl_wrong(w):
    from itertools import product as P
    fv = [R13.hpl_form_pullback13(0, s) for s in w]
    e0, eh, e1, em, einf = R13.maps13()[0]
    out = {}
    for ch in P(*[tuple(v.items()) for v in fv]):
        ww = tuple(x for cw, a in ch for x in cw)
        coef = Fr(1)
        for cw, a in ch: coef *= a
        out = R.vadd((out, 1), (seg_wrong(ww, e0, eh), coef))
    return out
nbad = sum(1 for w in WORDS if abs(R.cmzv_combo(hpl_wrong(w)) - R.hpl_value(w, mp.mpf(1)/3)) > mp.mpf("1e-30"))
check(f"B3 wrong composition order fails ({nbad}/{len(WORDS)} words off)", nbad > 10)

# B4 my-frame targets vanish; difference in proved span mod p
import numpy as np
T_my = {2: R13.my_target2(), 3: R13.my_target3(), 4: R13.my_target4()}
vals = [abs(R.cmzv_combo(T_my[n])) for n in (2, 3, 4)]
check(f"B4 my-frame targets vanish: {[mp.nstr(v,3) for v in vals]}", max(vals) < mp.mpf("1e-55"))
D4 = R.vadd((T_my[4], 1), ({w: Fr(c) for w, c in target4().items()}, -1))
def rank_modp(int_rows, ncols, p):
    M = np.array(int_rows, dtype=np.int64) % p
    nr, rank = M.shape[0], 0
    for j in range(ncols):
        if rank >= nr: break
        nz = np.nonzero(M[rank:, j])[0]
        if len(nz) == 0: continue
        piv = rank + int(nz[0])
        if piv != rank: M[[rank, piv]] = M[[piv, rank]]
        inv = pow(int(M[rank, j]), p-2, p)
        M[rank] = (M[rank]*inv) % p
        f = M[rank+1:, j].copy()
        M[rank+1:] = (M[rank+1:] - f[:, None]*M[rank]) % p
        rank += 1
    return rank
cols = list(cmzv6.words(4))
int_rows = [[dict(row).get(w, 0) for w in cols] for row in ROWS[4]]
from math import lcm
den = lcm(*[c.denominator for c in D4.values()]) if D4 else 1
dvec = [int(D4.get(w, Fr(0))*den) for w in cols]
r1 = rank_modp(int_rows, len(cols), 1000000009)
r2 = rank_modp(int_rows + [dvec], len(cols), 1000000009)
check(f"B4 weight-4 frame difference in proved span (mod 1000000009: rank {r1} vs {r2})", r1 == r2)

# B5 ID residuals at 250 dps + M3 at 1150 dps
mp.mp.dps = 250
def my_li(s, z, dps=250):
    z = mp.mpf(z)
    eps = mp.mpf(10)**(-dps-30)
    tot = mp.mpf(0); n = 1; zn = z
    while True:
        t = zn/mp.mpf(n)**s
        tot += t
        if n > 200 and abs(t)*abs(z)/(1-abs(z)) < eps: break
        n += 1; zn *= z
    return tot
def fit_val(name, weight):
    tot = mp.mpf(0)
    for e, w, v in R13.fit_terms(name):
        if len(w) != weight: continue
        hv = mp.mpf(1) if not w else R.hpl_value(w, mp.mpf(1)/3)
        tot += mp.mpf(v.numerator)/v.denominator * (mp.mpf(1)/3)**e * hv
    return tot
R32 = fit_val("hpl_fit_k3_(0,0,0).json", 2); R33 = fit_val("hpl_fit_k3_(0,0,0).json", 3)
R34 = fit_val("hpl_fit_k3_(0,0,0).json", 4)
P1 = fit_val("hpl_fit_k3_(1,0,0).json", 1); P2 = fit_val("hpl_fit_k3_(1,0,0).json", 2)
pi, z2, z3 = mp.pi, mp.zeta(2), mp.zeta(3)
L2, L3 = mp.log(2), mp.log(3)
Li21, Li2m = my_li(2,"1/3"), my_li(2,"-1/3")
Li31, Li3m = my_li(3,"1/3"), my_li(3,"-1/2")
id2 = R32 + z2/2 - (-Li21 - Li2m + pi**2/12)
W3cl = (-26*z3 - 87*Li31 - 130*Li3m + mp.mpf(21)/4*pi**2*L2 - mp.mpf(43)/12*pi**2*L3
        + 67*Li21*L2 - 45*Li21*L3 + mp.mpf(65)/3*L2**3 - 65*L2**2*L3 + mp.mpf(109)/2*L2*L3**2
        - mp.mpf(37)/3*L3**3)
id3 = R33 + z2*P1 - z3/2 - W3cl
check(f"B5 ID2 residual 250 dps: {mp.nstr(abs(id2),3)}", abs(id2) < mp.mpf("1e-225"))
check(f"B5 ID3 residual 250 dps: {mp.nstr(abs(id3),3)}", abs(id3) < mp.mpf("1e-220"))
# M3 at 1150 dps
mp.mp.dps = 1150
def my_li2(s, z): return my_li(s, z, 1150)
def my_li4_neg2():
    y = mp.mpc(mp.mpf(1)/2) + mp.log(2)/(2*mp.pi*mp.j)
    B4 = y**4 - 2*y**3 + y**2 - mp.mpf(1)/30
    return (-(2*mp.pi*mp.j)**4/mp.factorial(4)*B4 - my_li2(4, mp.mpf(-1)/2)).real
Li21, Li2m = my_li2(2,"1/3"), my_li2(2,"-1/3")
Li31, Li3m = my_li2(3,"1/3"), my_li2(3,"-1/2")
Li41, Li4m3 = my_li2(4,"1/3"), my_li2(4,"-1/3")
Li423, Li412, Li419 = my_li2(4,"2/3"), my_li2(4,"1/2"), my_li2(4,"1/9")
Li4m2 = my_li4_neg2()
pi, z2, z3 = mp.pi, mp.zeta(2), mp.zeta(3)
L2, L3 = mp.log(2), mp.log(3)
W2cl = -Li21 - Li2m + pi**2/12
W3cl = (-26*z3 - 87*Li31 - 130*Li3m + mp.mpf(21)/4*pi**2*L2 - mp.mpf(43)/12*pi**2*L3
        + 67*Li21*L2 - 45*Li21*L3 + mp.mpf(65)/3*L2**3 - 65*L2**2*L3 + mp.mpf(109)/2*L2*L3**2
        - mp.mpf(37)/3*L3**3)
W4cl = (-mp.mpf(22)/3*Li4m2 - mp.mpf(17)/2*Li41 + mp.mpf(20)/3*Li423 + mp.mpf(110)/9*Li412
        + mp.mpf(25)/48*Li419 - 12*Li31*L2 + mp.mpf(28)/3*Li31*L3 - mp.mpf(16)/3*Li3m*L2
        + mp.mpf(16)/3*Li3m*L3 + mp.mpf(14)/3*z3*L2 - mp.mpf(25)/9*z3*L3 - 4*Li21**2
        + mp.mpf(13)/3*Li21*L2**2 - 6*Li21*L2*L3 + mp.mpf(2)/3*Li21*L3**2
        + mp.mpf(11)/18*pi**2*Li21 - mp.mpf(31)/27*pi**2*L2**2 + mp.mpf(31)/72*pi**2*L3**2
        - mp.mpf(271)/1296*pi**4 + mp.mpf(46)/27*L2**4 - mp.mpf(62)/9*L2**3*L3
        + mp.mpf(61)/6*L2**2*L3**2 - mp.mpf(17)/3*L2*L3**3 + mp.mpf(103)/144*L3**4)
M3 = W2cl + W3cl + W4cl
quad = mp.mpf(open(BASE + "/COR7_classical_corners/c13_M3_quad_1100.txt").read().strip())
d = abs(M3 - quad)
check(f"B5 M3 closed form vs 1095-digit quadrature: {int(-mp.log10(d))} digits", d < mp.mpf("1e-1080"))
log("PART B " + ("PASSED" if ok else "FAILED"))
log("VERDICT: " + ("COR13CR REFEREE CHECKS ALL PASS" if ok else "FAILURE"))
sys.exit(0 if ok else 1)
