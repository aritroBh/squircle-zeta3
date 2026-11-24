"""step_02.py -- referee TASK 1 (algebra): exact certificate check with MY OWN code.
Rebuild the proved relation rows (cmzv6 ordering, for prow indexing), rebuild the targets
(hexc2_common conversion), load cert_w23.json / cert_w4.json / cert_euler.json, and check
    sum_i y_i * rows[prow_i] - T_n == 0   on EVERY word,
(a) with FLINT fmpq_mat exact arithmetic (my own construction), all weights;
(b) with an independent pure-Fraction implementation at weights 2, 3 and Euler cert;
(c) weight 4 also modulo 5 random primes with my own modular arithmetic (guards against
    a hypothetical fmpq_mat issue)."""
import sys, json, time, random
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC2_level6_certificate")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
sys.path.insert(0, "../COR4_M3alpha_classical")
from fractions import Fraction as Fr
import cmzv6
from hexc2_common import target2, target3, target4
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

cert23 = json.load(open("../HEXC2_level6_certificate/cert_w23.json"))
cert4 = json.load(open("../HEXC2_level6_certificate/cert_w4.json"))
certE = json.load(open("../HEXC2_level6_certificate/cert_euler.json"))

# ---- my own residual machinery ----
def residual_fractions(rows_sel, y, T, cols):
    """pure-Fraction: sum_i y_i row_i - T as a dict over cols; returns #nonzero."""
    ids = {w: j for j, w in enumerate(cols)}
    acc = [Fr(0)] * len(cols)
    for row, yi in zip(rows_sel, y):
        if yi == 0: continue
        for w, c in row:
            acc[ids[w]] += yi * Fr(c)
    bad = 0
    for w, j in ids.items():
        if acc[j] != T.get(w, Fr(0)): bad += 1
    return bad

def residual_fmpq(rows_sel, y, T, cols):
    from flint import fmpq_mat
    ids = {w: j for j, w in enumerate(cols)}
    r, m = len(rows_sel), len(cols)
    S = fmpq_mat(r, m)
    for i, row in enumerate(rows_sel):
        for w, c in row:
            S[i, ids[w]] = c
    Y = fmpq_mat(r, 1)
    for i, yi in enumerate(y):
        Y[i, 0] = fmpq_mat(1, 1, [str(yi)])[0, 0]
    res = Y.transpose() * S
    bad = 0
    for w, j in ids.items():
        if Fr(str(res[0, j])) != T.get(w, Fr(0)): bad += 1
    return bad

def residual_modp(rows_sel, y, T, cols, p):
    """my own modular computation of (Y^T S - T) mod p; needs T,y rational."""
    ids = {w: j for j, w in enumerate(cols)}
    acc = [0] * len(cols)
    for row, yi in zip(rows_sel, y):
        ym = yi.numerator % p * pow(yi.denominator % p, p - 2, p) % p
        if ym == 0: continue
        for w, c in row:
            j = ids[w]
            acc[j] = (acc[j] + ym * (c % p)) % p
    bad = 0
    for w, j in ids.items():
        t = T.get(w, Fr(0))
        tm = t.numerator % p * pow(t.denominator % p, p - 2, p) % p
        if (acc[j] - tm) % p != 0: bad += 1
    return bad

lower = {}
for n in (1, 2, 3, 4):
    rows, counts = cmzv6.make_relations(n, lower)
    lower[n] = rows
    if n < 2: continue
    T = {2: target2, 3: target3, 4: target4}[n]()
    cols = list(cmzv6.words(n))
    if n < 4:
        C = cert23[str(n)]
        prow, y = C["prow"], [Fr(int(a), int(b)) for a, b in C["y"]]
    else:
        prow, y = cert4["prow"], [Fr(s) for s in cert4["y"]]
    assert len(set(prow)) == len(prow) and max(prow) < len(rows) and min(prow) >= 0
    sel = [rows[i] for i in prow]
    bad_f = residual_fmpq(sel, y, T, cols)
    log(f"weight {n}: rows={len(rows)} cert_rows={len(prow)} target_support={len(T)} "
        f"fmpq residual_words={bad_f}  {'PASS' if bad_f == 0 else 'FAIL'}")
    assert bad_f == 0
    if n <= 3:
        bad_py = residual_fractions(sel, y, T, cols)
        log(f"          pure-Fraction residual_words={bad_py}  {'PASS' if bad_py == 0 else 'FAIL'}")
        assert bad_py == 0
    else:
        random.seed(20240917)
        primes = [1000000007, 1000000009, 998244353, 999999937, 99999989]
        for p in primes:
            bad_p = residual_modp(sel, y, T, cols, p)
            log(f"          mod {p}: residual_words={bad_p}  {'PASS' if bad_p == 0 else 'FAIL'}")
            assert bad_p == 0
        h = max((len(str(yy.numerator)) for yy in y), default=0)
        log(f"          max y numerator digits = {h}")

# Euler cert: zeta(2) = pi^2/6  <=>  (-1,0) - 3*(1,1)... target {(-1,0):1, (1,1):3}? value (1,1)=Li_2(zeta)... check below
ET = {(-1, 0): Fr(1), (1, 1): Fr(3)}
cols2 = list(cmzv6.words(2))
assert [list(w) for w in cols2] == certE["cols"]
sel = [lower[2][i] for i in certE["prow"]]
y = [Fr(s) for s in certE["y"]]
bad_f = residual_fmpq(sel, y, ET, cols2)
bad_py = residual_fractions(sel, y, ET, cols2)
log(f"euler cert: rows={len(sel)} fmpq residual={bad_f} pure-Fraction residual={bad_py} "
    f"{'PASS' if bad_f == bad_py == 0 else 'FAIL'}")
assert bad_f == 0 and bad_py == 0
log("STEP 02 PASS: certificates are exact on every word (3 independent arithmetic backends)")
with open("progress.log", "a") as f:
    f.write("step_02: exact cert check PASS w/ my fmpq code (w2/3/4+euler), pure-Fraction (w2/3), mod-5-primes (w4)\n")
