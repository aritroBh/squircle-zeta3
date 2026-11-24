import time, sys, json, os, random
sys.dont_write_bytecode = True
base = '<project>/research/08_theory'
sys.path.insert(0, base + '/HEXC1_hexagon_certificate')
sys.path.insert(0, base + '/COR13C_third_corner_certificate')
sys.path.insert(0, base + '/COR7_classical_corners')
from fractions import Fraction as Fr
import cmzv6
from cor13c_common import target2, target3, target4
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

D = base + '/COR13C_third_corner_certificate'
cert23 = json.load(open(D + '/cert13_w23.json'))
cert4 = json.load(open(D + '/cert13_w4.json'))
certE = json.load(open(D + '/cert13_euler.json'))

def residual_fractions(rows_sel, y, T, cols):
    ids = {w: j for j, w in enumerate(cols)}
    acc = [Fr(0)] * len(cols)
    for row, yi in zip(rows_sel, y):
        if yi == 0: continue
        for w, c in row:
            acc[ids[w]] += yi * Fr(c)
    return sum(1 for w, j in ids.items() if acc[j] != T.get(w, Fr(0)))

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
    return sum(1 for w, j in ids.items() if Fr(str(res[0, j])) != T.get(w, Fr(0)))

def residual_modp(rows_sel, y, T, cols, p):
    ids = {w: j for j, w in enumerate(cols)}
    acc = [0] * len(cols)
    for row, yi in zip(rows_sel, y):
        ym = yi.numerator % p * pow(yi.denominator % p, p - 2, p) % p
        if ym == 0: continue
        for w, c in row:
            acc[ids[w]] = (acc[ids[w]] + ym * (c % p)) % p
    bad = 0
    for w, j in ids.items():
        t = T.get(w, Fr(0))
        tm = t.numerator % p * pow(t.denominator % p, p - 2, p) % p
        if (acc[j] - tm) % p: bad += 1
    return bad

lower = {}
supports = {}
for n in (1, 2, 3, 4):
    rows, counts = cmzv6.make_relations(n, lower)
    lower[n] = rows
    if n < 2: continue
    T = {2: target2, 3: target3, 4: target4}[n]()
    supports[n] = len(T)
    cols = list(cmzv6.words(n))
    if n < 4:
        C = cert23[str(n)]
        prow, y = C["prow"], [Fr(int(a), int(b)) for a, b in C["y"]]
        assert [list(w) for w in cols] == C["cols"]
    else:
        prow, y = cert4["prow"], [Fr(s) for s in cert4["y"]]
        assert [list(w) for w in cols] == cert4["cols"]
    assert len(set(prow)) == len(prow) and max(prow) < len(rows) and min(prow) >= 0
    sel = [rows[i] for i in prow]
    # row supports must be inside admissible words
    for row in sel:
        for w, c in row:
            assert w in set(cols) or (w[0] != 0 and w[-1] != -1)
    bad_f = residual_fmpq(sel, y, T, cols)
    log(f"weight {n}: rows={len(rows)} cert_rows={len(prow)} target_support={len(T)} fmpq residual_words={bad_f}")
    assert bad_f == 0
    if n <= 3:
        bad_py = residual_fractions(sel, y, T, cols)
        log(f"          pure-Fraction residual_words={bad_py}")
        assert bad_py == 0
    else:
        for p in (1000000007, 1000000009, 998244353, 999999937, 99999989):
            bad_p = residual_modp(sel, y, T, cols, p)
            assert bad_p == 0, (p, bad_p)
        log(f"          mod 5 primes residual_words=0; max y numerator digits = {max(len(str(yy.numerator)) for yy in y)}")
        log(f"          T4 coeff of (1,1,1,1) = {T.get((1,1,1,1), Fr(0))}")

# Euler cert
ET = {(-1, 0): Fr(1), (1, 1): Fr(3)}
cols2 = list(cmzv6.words(2))
sel = [lower[2][i] for i in certE["prow"]]
y = [Fr(int(a), int(b)) for a, b in certE["y"]]
bad_f = residual_fmpq(sel, y, ET, cols2); bad_py = residual_fractions(sel, y, ET, cols2)
log(f"euler cert: {len(sel)} rows; fmpq residual={bad_f}, pure-Fraction residual={bad_py}")
assert bad_f == 0 and bad_py == 0
log(f"target supports (mine from cor13c_common): {supports}")
log("STEP PASS: all certificates exact on every word, three arithmetic backends")
with open('progress.log','a') as f:
    f.write("referee step(cert exact): w2/w3/w4+euler residual 0 on all words; fmpq + pure-Fraction (w2/3,euler) + mod 5 primes (w4)\n")
