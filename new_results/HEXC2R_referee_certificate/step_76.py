
verifier = '''"""verify_hexc2r.py -- HEXC2R referee's independent verifier. Runtime target < 5 min.

PART A (exact): certificate algebra with MY OWN fmpq code: sum y_i rows[prow_i] - T == 0 on
  every word, weights 2/3/4 + Euler cert + the referee's Psi2 bonus cert.
PART B (numeric, >= 60 dps): >= 50 relation rows across all families vanish under MY OWN
  CMZV evaluator (truncated iterated sums; cross-checked against HEXC2's evaluators).
PART C (conversion): my conversion EXACTLY equals hexc2_common on all 38 HPL words + 9 Li
  atoms; 10-word numeric audit at 65 dps; 36 pulled-back 1-form identities EXACT over Q(zeta).
PART D (numbers): M3(2pi/3) from the closed form with my own Li2/Li3/Li4 (+ Jonquiere for
  Li4(-2)) at 260 dps vs COR5's value; quotient dims weights 1-3 mod p; Psi2 closed form.
Full audits (all 1952 cert rows numerically, weight-4 rank, 1094-digit M3) are in
steps/step_04.py, steps/step_09.py, steps/step_08.py (offline evidence).
"""
import os, sys, json, re, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, HERE + "/../HEXC2_level6_certificate")
sys.path.insert(0, HERE + "/../HEXC1_hexagon_certificate")
sys.path.insert(0, HERE + "/../COR4_M3alpha_classical")
from fractions import Fraction as Fr
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)
ok = True
def require(cond, msg):
    global ok
    ok = ok and cond
    log(("  [PASS] " if cond else "  [FAIL] ") + msg)
    assert cond, msg

# ================= PART A (exact) =================
log("PART A: exact certificate verification (my fmpq code)")
import cmzv6
from hexc2_common import target2, target3, target4
from flint import fmpq_mat, fmpq

cert23 = json.load(open(HERE + "/../HEXC2_level6_certificate/cert_w23.json"))
cert4 = json.load(open(HERE + "/../HEXC2_level6_certificate/cert_w4.json"))
certE = json.load(open(HERE + "/../HEXC2_level6_certificate/cert_euler.json"))

def exact_check(rows, cols, prow, y, T):
    ids = {w: i for i, w in enumerate(cols)}
    sel = [rows[i] for i in prow]
    S = fmpq_mat(len(sel), len(cols))
    for i, row in enumerate(sel):
        for w, c in row:
            S[i, ids[w]] = c
    Y = fmpq_mat(len(sel), 1)
    for i, c in enumerate(y):
        Y[i, 0] = fmpq(str(c))
    res = Y.transpose() * S
    return sum(1 for w, j in ids.items() if Fr(str(res[0, j])) != T.get(w, Fr(0)))

lower = {}
for n in (1, 2, 3, 4):
    rows, counts = cmzv6.make_relations(n, lower)
    lower[n] = rows
    if n < 2: continue
    T = {2: target2, 3: target3, 4: target4}[n]()
    cols = list(cmzv6.words(n))
    if n < 4:
        C = cert23[str(n)]; prow, y = C["prow"], [Fr(int(a), int(b)) for a, b in C["y"]]
    else:
        prow, y = cert4["prow"], [Fr(s) for s in cert4["y"]]
    bad = exact_check(rows, cols, prow, y, T)
    require(bad == 0, f"weight {n}: exact cert residual 0 on {len(cols)} words ({len(prow)} rows)")
ET = {(-1, 0): Fr(1), (1, 1): Fr(3)}
bad = exact_check(lower[2], list(cmzv6.words(2)), certE["prow"], [Fr(s) for s in certE["y"]], ET)
require(bad == 0, "Euler cert: zeta(2) = pi^2/6 exact in weight-2 span")
# referee's Psi2 bonus cert
import reflib as R
CP = json.load(open(HERE + "/cert_psi2.json"))
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
def _prod(*vs):
    out = {(): Fr(1)}
    for v in vs: out = R.mul(out, v)
    return out
hpl = {}
for key, val in json.load(open(G1 + "/hpl_fit_k3_(1,0,0).json")):
    m = re.match(r"c\\^(-?\\d+) H\\[(.*)\\]", key)
    if not m: continue
    e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
    if len(w) != 2: continue
    hpl = R.vadd((hpl, 1), (R.my_hpl_to_cmzv(w), Fr(val) * Fr(1, 2) ** e))
PI2 = {w: -c for w, c in R.my_atom_pi2().items()}
closed = R.vadd((PI2, Fr(3, 32)), (R.my_li_atom(2, Fr(1, 3)), Fr(-9, 8)),
                (_prod(R.my_atom_log2(), R.my_atom_log3()), Fr(3, 4)),
                (_prod(R.my_atom_log3(), R.my_atom_log3()), Fr(-9, 16)))
bad = exact_check(lower[2], [tuple(w) for w in CP["cols"]], CP["prow"], [Fr(s) for s in CP["y"]],
                  R.vadd((hpl, 1), (closed, -1)))
require(bad == 0, "referee bonus: Psi2(1/2) closed form exact (28 rows, residual 0)")

# ================= PART B (numeric, 62 dps) =================
log("PART B: relation-row numeric sample with MY evaluator (62 dps)")
import mpmath as mp
mp.mp.dps = 62
sys.path.insert(0, HERE + "/../COR4_M3alpha_classical")
import numeric_cmzv as NC
NC.set_dps(60)
worst3 = 0
for w in [(-1, 0), (1, 3, -1, 2), (2, 0, 4, 1)]:
    worst3 = max(worst3, abs(R.cmzv_value(w) - NC.cmzv_value(w)))
require(worst3 < mp.mpf("1e-50"), f"my evaluator vs HEXC2 numeric_cmzv (3 words): {mp.nstr(worst3, 3)}")

def tag_families(n, lower):
    fam = {}
    def emit(row, f):
        r = cmzv6.primitive(row)
        if r and r not in fam: fam[r] = f
    if n >= 2:
        for a in range(1, n):
            b = n - a
            for u in cmzv6.words(a, False):
                u_adm = u[0] != 0
                if u_adm and a > b: continue
                for v in cmzv6.words(b):
                    if u_adm and a == b and u > v: continue
                    emit(cmzv6.ds(u, v), 'DS')
    for d in (2, 3, 6):
        for row in cmzv6.distribution(n, d): emit(row, f'dist{d}')
    for k, basis in lower.items():
        if k >= n: continue
        for row in basis:
            for w in cmzv6.words(n - k):
                emit(cmzv6.mul(dict(row), {w: 1}), f'lift{k}')
    return fam
import random
random.seed(20250917)
tested = 0; overall = 0; famseen = set()
for n in (2, 3, 4):
    rows = lower[n]
    fam = tag_families(n, lower)
    prow = cert23[str(n)]["prow"] if n < 4 else cert4["prow"]
    # all cert rows at weight 2; stratified sample at weights 3,4 covering every family
    if n == 2:
        idx = list(prow)
    else:
        idx = []
        per = {}
        cand = list(prow); random.shuffle(cand)
        for i in cand:
            f = fam[rows[i]]
            if per.get(f, 0) < (5 if n == 3 else 4):
                per[f] = per.get(f, 0) + 1; idx.append(i)
    worst = 0
    for i in idx:
        val = mp.mpc(0)
        for w, c in rows[i]:
            val += int(c) * R.cmzv_value(w)
        worst = max(worst, abs(val)); famseen.add(fam[rows[i]])
    tested += len(idx); overall = max(overall, worst)
    log(f"  weight {n}: {len(idx)} rows, worst {mp.nstr(worst, 3)}")
require(tested >= 50 and overall < mp.mpf("1e-55"),
        f"{tested} rows ({sorted(famseen)}) vanish at 62 dps, worst {mp.nstr(overall, 3)}")

# ================= PART C (conversion) =================
log("PART C: conversion checks")
import hexc2_common as HC
def load_fit(name):
    D = json.load(open(G1 + "/" + name)); out = []
    for key, val in D:
        m = re.match(r"c\\^(-?\\d+) H\\[(.*)\\]", key)
        if key == "1": e, w = 0, ()
        elif key == "c": e, w = 1, ()
        else:
            e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
        out.append((e, w))
    return out
WORDS = sorted({w for name in ["hpl_fit_k3_(0,0,0).json", "hpl_fit_k3_(1,0,0).json"]
                for e, w in load_fit(name) if w}, key=lambda x: (len(x), x))
same = all(R.my_hpl_to_cmzv(w) == HC.hpl_to_cmzv(w) for w in WORDS)
require(same, f"my conversion == hexc2 EXACTLY on all {len(WORDS)} HPL words")
LI = [(2, "1/3"), (3, "1/3"), (3, "-1/2"), (4, "1/3"), (4, "-1/3"), (4, "2/3"),
      (4, "1/2"), (4, "-1/2"), (4, "-2")]
same = True
for n, z in LI:
    frame = tuple(None if p is HC.INF else p for p in HC.cross_table()[Fr(z)][0])
    same &= (R.my_li_atom(n, Fr(z), frame=frame) == HC.li_atom(n, Fr(z)))
require(same, f"my Li-atom conversion == hexc2 EXACTLY on all {len(LI)} atoms")
mp.mp.dps = 65
sample = WORDS[:3] + WORDS[13:15] + WORDS[24:26] + WORDS[-3:]
worst = 0
for w in sample:
    worst = max(worst, abs(R.cmzv_combo(R.my_hpl_to_cmzv(w)) - R.hpl_value(w, mp.mpf(1)/2)))
require(worst < mp.mpf("1e-55"), f"{len(sample)}-word numeric audit of conversions: {mp.nstr(worst, 3)}")
worst = 0
for n, z in LI[:4]:
    worst = max(worst, abs(R.cmzv_combo(R.my_li_atom(n, Fr(z))) - mp.polylog(n, mp.mpf(z))))
require(worst < mp.mpf("1e-55"), f"4 Li atoms numeric: {mp.nstr(worst, 3)}")
# exact 1-form pullbacks (my Q(zeta) arithmetic)
def qadd(x, y): return (x[0] + y[0], x[1] + y[1])
def qmul(x, y):
    a, b = x; c, d = y
    return (a*c - b*d, a*d + b*c + b*d)
def qdiv(x, y):
    a, b = y; nn = a*a + a*b + b*b
    return qmul(x, (Fr(a + b, nn), Fr(-b, nn)))
def qpow(x, n):
    out = (Fr(1), Fr(0))
    for _ in range(n): out = qmul(out, x)
    return out
def qsub(x, y): return (x[0] - y[0], x[1] - y[1])
ZETA = (Fr(0), Fr(1)); MU6 = [qpow(ZETA, k) for k in range(6)]; Z0 = (Fr(0), Fr(0))
def pm(P, Q):
    out = [Z0] * (len(P) + len(Q) - 1)
    for i, a in enumerate(P):
        for j, b in enumerate(Q): out[i + j] = qadd(out[i + j], qmul(a, b))
    return out
def pa(P, Q):
    n = max(len(P), len(Q))
    return [qadd(P[i] if i < len(P) else Z0, Q[i] if i < len(Q) else Z0) for i in range(n)]
def ps(P, c): return [qmul(a, c) for a in P]
def trim(P):
    while len(P) > 1 and P[-1] == Z0: P = P[:-1]
    return P
def form_to_poly(lettervec):
    terms = [(c, w[0]) for w, c in lettervec.items()]
    den = [(Fr(1), Fr(0))]
    if any(j == -1 for c, j in terms): den = pm(den, [Z0, (Fr(1), Fr(0))])
    for c, j in terms:
        if j >= 0: den = pm(den, [(Fr(1), Fr(0)), (-MU6[j][0], -MU6[j][1])])
    num = [Z0]
    for c, j in terms:
        p = [(Fr(1), Fr(0))]
        for c2, j2 in terms:
            if j2 == j: continue
            p = pm(p, [Z0, (Fr(1), Fr(0))] if j2 == -1 else [(Fr(1), Fr(0)), (-MU6[j2][0], -MU6[j2][1])])
        cc = MU6[j] if j >= 0 else (Fr(1), Fr(0))
        num = pa(num, ps(p, qmul((Fr(c.numerator, c.denominator), Fr(0)), cc)))
    return num, den
def mob3(f0, f1, fi):
    A = fi; C = (Fr(1), Fr(0))
    D = qdiv(qsub(fi, f1), qsub(f1, f0))
    return (A, qmul(f0, D), C, D)
def invm(M):
    A, B, C, D = M
    return (D, (-B[0], -B[1]), (-C[0], -C[1]), A)
cnt = 0
for mi, E in enumerate(R.MAPS6):
    A, B, C, D = invm(mob3(MU6[E[0]], MU6[E[2]], MU6[E[4]]))
    for sig in (0, 1, -1):
        if sig == 0:
            Rn = [qsub(qmul(A, D), qmul(B, C))]; Rd = pm([B, A], [D, C])
        else:
            A2, B2, sg = (qsub(A, C), qsub(B, D), (Fr(-1), Fr(0))) if sig == 1 else (qadd(A, C), qadd(B, D), (Fr(1), Fr(0)))
            Rn = [qmul(sg, qsub(qmul(A2, D), qmul(B2, C)))]; Rd = pm([B2, A2], [D, C])
        Cn, Cd = form_to_poly(R.hpl_form_pullback(mi, sig))
        assert trim(pm(Rn, Cd)) == trim(pm(Cn, Rd)), (mi, sig)
        cnt += 1
require(cnt == 36, "36 pulled-back 1-form identities EXACT over Q(zeta)")

# ================= PART D (numbers) =================
log("PART D: M3(2pi/3) from the closed form (260 dps, own polylogs)")
mp.mp.dps = 260
def my_pl(s, z):
    z = mp.mpc(z); az = abs(z); eps = mp.mpf(10) ** (-290)
    tot = mp.mpc(0); n = 1; zn = z
    while True:
        t = zn / mp.mpf(n) ** s
        tot += t
        if n > 100 and abs(t) * az / (1 - az) < eps: break
        n += 1; zn *= z
    return tot
def my_li4_neg2():
    y = mp.mpc(mp.mpf(1)/2) + mp.log(2) / (2 * mp.pi * mp.j)
    B4 = y**4 - 2*y**3 + y**2 - mp.mpf(1)/30
    return -(2 * mp.pi * mp.j) ** 4 / 24 * B4 - my_pl(4, mp.mpf(-1)/2)
require(abs(my_li4_neg2() - mp.polylog(4, -2)) < mp.mpf("1e-240"), "Li4(-2) Jonquiere vs mpmath")
pi, z2, z3, L2, L3 = mp.pi, mp.zeta(2), mp.zeta(3), mp.log(2), mp.log(3)
Li21 = my_pl(2, mp.mpf(1)/3); Li31 = my_pl(3, mp.mpf(1)/3); Li3m = my_pl(3, mp.mpf(-1)/2)
Li41 = my_pl(4, mp.mpf(1)/3); Li4m3 = my_pl(4, mp.mpf(-1)/3); Li423 = my_pl(4, mp.mpf(2)/3)
Li412 = my_pl(4, mp.mpf(1)/2); Li419 = my_pl(4, mp.mpf(1)/9); Li4m2 = my_li4_neg2()
W2 = L2**2/2 + Li21 + mp.log(mp.mpf(3)/2)**2/2
W3 = mp.mpf(-15)/4*z3 - 3*Li31 - 8*Li3m - pi**2*L2/12 + 4*Li21*L2 - mp.mpf(3)/2*Li21*L3 \\
     + L2**3/3 - mp.mpf(5)/2*L2**2*L3 + 2*L2*L3**2 - L3**3/4
F4 = (-21*Li4m2 - mp.mpf(65)/8*Li41 + mp.mpf(13)/2*Li423 + mp.mpf(137)/6*Li412 + mp.mpf(39)/64*Li419
      + mp.mpf(15)/2*Li31*L2 + mp.mpf(37)/2*Li3m*L2 - mp.mpf(5)/4*Li3m*L3 + mp.mpf(119)/16*z3*L2
      - mp.mpf(37)/32*z3*L3 + Li21**2 - mp.mpf(23)/4*Li21*L2**2 + mp.mpf(7)/4*Li21*L2*L3
      + Li21*L3**2 - pi**2*Li21/6 - mp.mpf(455)/144*pi**2*L2**2 + mp.mpf(9)/8*pi**2*L2*L3
      - mp.mpf(7)/32*pi**2*L3**2 - mp.mpf(949)/2160*pi**4 - mp.mpf(211)/144*L2**4
      + mp.mpf(101)/24*L2**3*L3 - mp.mpf(15)/8*L2**2*L3**2 - mp.mpf(35)/24*L2*L3**3
      + mp.mpf(113)/192*L3**4)
psi2 = mp.mpf(3)/32*pi**2 - mp.mpf(9)/8*Li21 + mp.mpf(3)/4*L2*L3 - mp.mpf(9)/16*L3**2
W4 = F4 + z2 * psi2 + mp.mpf(3)/8 * z3 * L3
M3 = W4 + W3 + W2
m3ref = open(HERE + "/../COR5_corner_classicality/M3_2pi3_1100.txt").read().strip()[:262]
d = abs(M3 - mp.mpf(m3ref))
require(d < mp.mpf("1e-250"), f"M3 closed form vs COR5 value: |diff| = {mp.nstr(d, 3)}")
# quotient dims weights 1-3 mod p (weight 4 in steps/step_09.py)
import numpy as np
def rank_p(int_rows, ncols, p):
    M = np.array(int_rows, dtype=np.int64) % p
    nr, rank = M.shape[0], 0
    for j in range(ncols):
        if rank >= nr: break
        nz = np.nonzero(M[rank:, j])[0]
        if len(nz) == 0: continue
        piv = rank + int(nz[0])
        M[[rank, piv]] = M[[piv, rank]]
        M[rank] = (M[rank] * pow(int(M[rank, j]), p - 2, p)) % p
        f = M[rank + 1:, j].copy()
        M[rank + 1:] = (M[rank + 1:] - f[:, None] * M[rank]) % p
        rank += 1
    return rank
quot = {}
for n in (1, 2, 3):
    rows = lower[n]; cols = list(cmzv6.words(n)); ids = {w: j for j, w in enumerate(cols)}
    ir = [[0] * len(cols) for _ in range(len(rows))]
    for i, row in enumerate(rows):
        for w, c in row: ir[i][ids[w]] = c
    quot[n] = len(cols) - rank_p(ir, len(cols), 1000000007)
require(quot == {1: 3, 2: 8, 3: 23}, f"quotient dims weights 1-3 mod p: {quot} (w4 = 69 in step_09)")
log("VERDICT: " + ("HEXC2R: ALL REFEREE CHECKS PASS — certificate VALID" if ok else "FAILURES PRESENT"))
assert ok
'''
with open("verify_hexc2r.py", "w") as f:
    f.write(verifier)
import subprocess, time
t0 = time.monotonic()
r = subprocess.run(["nice", "-n", "10", "<project>/.venv/bin/python", "verify_hexc2r.py"], capture_output=True, text=True, timeout=1500)
print(f"elapsed: {time.monotonic()-t0:.1f}s")
print(r.stdout)
print(r.stderr[-2500:] if r.returncode else "OK")
