"""step_11.py -- referee BONUS certificate: Psi^{(2)}(1/2) in closed form, proved exactly.
From the G1 fit (weight-2 terms) and the PSLQ guess
  Psi2(1/2) = (3/32) pi^2 - (9/8) Li_2(1/3) + (3/4) log2 log3 - (9/16) log^2 3,
build the level-6 weight-2 target with MY conversion and solve exactly for multipliers on
the 28-row weight-2 basis (cert rows), verify on all 36 words exactly, and numerically."""
import sys, json, re, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC2_level6_certificate")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
import mpmath as mp
from fractions import Fraction as Fr
import reflib as R
import cmzv6
from flint import fmpq_mat, fmpq
import numpy as np
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
def prod(*vecs):
    out = {(): Fr(1)}
    for v in vecs: out = R.mul(out, v)
    return out

# Psi2 HPL side (my conversion)
hpl = {}
for key, val in json.load(open(G1 + "/hpl_fit_k3_(1,0,0).json")):
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", key)
    if not m: continue
    e = int(m.group(1))
    w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
    if len(w) != 2: continue
    hpl = R.vadd((hpl, 1), (R.my_hpl_to_cmzv(w), Fr(val) * Fr(1, 2) ** e))
# closed form side (my atoms)
PI2 = {w: -c for w, c in R.my_atom_pi2().items()}     # pi^2 = -(i pi)^2
closed = R.vadd(
    (PI2, Fr(3, 32)),
    (R.my_li_atom(2, Fr(1, 3)), Fr(-9, 8)),
    (prod(R.my_atom_log2(), R.my_atom_log3()), Fr(3, 4)),
    (prod(R.my_atom_log3(), R.my_atom_log3()), Fr(-9, 16)),
)
T = R.vadd((hpl, 1), (closed, -1))
log(f"target Psi2 - closed: support {len(T)}")

# numeric first
mp.mp.dps = 65
v = R.cmzv_combo(T)
log(f"|target value| = {mp.nstr(abs(v), 3)} (65 dps)")
assert abs(v) < mp.mpf("1e-55")

# exact solve on the weight-2 cert rows (rank 28)
lower = {}
rows1, _ = cmzv6.make_relations(1, lower); lower[1] = rows1
rows, _ = cmzv6.make_relations(2, lower)
cert = json.load(open("../HEXC2_level6_certificate/cert_w23.json"))["2"]
prow = cert["prow"]
sel = [rows[i] for i in prow]
cols = list(cmzv6.words(2))
ids = {w: j for j, w in enumerate(cols)}
# choose 28 pivot words mod p with my own elimination (record pivots)
p = 1000000007
M = np.zeros((len(sel), len(cols)), dtype=np.int64)
for i, row in enumerate(sel):
    for w, c in row:
        M[i, ids[w]] = c % p
nr, rank = M.shape[0], 0
pivots = []
for j in range(len(cols)):
    col = M[rank:, j]
    nz = np.nonzero(col)[0]
    if len(nz) == 0: continue
    piv = rank + int(nz[0])
    M[[rank, piv]] = M[[piv, rank]]
    inv = pow(int(M[rank, j]), p - 2, p)
    M[rank] = (M[rank] * inv) % p
    f = M[rank + 1:, j].copy()
    M[rank + 1:] = (M[rank + 1:] - f[:, None] * M[rank]) % p
    pivots.append(j)
    rank += 1
    if rank == nr: break
assert rank == len(sel) == 28, (rank, len(sel))
# exact solve: square subsystem on pivot columns
S = fmpq_mat(28, 28)
for i, row in enumerate(sel):
    for w, c in row:
        if ids[w] in pivots:
            S[pivots.index(ids[w]), i] = c   # S = rows^T restricted to pivot words: 28 x 28
b = fmpq_mat(28, 1)
for k, j in enumerate(pivots):
    b[k, 0] = fmpq(str(T.get(cols[j], Fr(0))))
X = S.solve(b)
y = [Fr(str(X[k, 0])) for k in range(28)]
# verify on ALL words exactly
acc = {}
for yi, row in zip(y, sel):
    for w, c in row:
        acc[w] = acc.get(w, Fr(0)) + yi * Fr(c)
bad = [w for w in cols if acc.get(w, Fr(0)) != T.get(w, Fr(0))]
log(f"exact solve on 28 pivot words; residual on all 36 words: {len(bad)}  "
    f"{'PASS' if not bad else 'FAIL'}")
assert not bad
h = max(len(str(yy.numerator)) for yy in y)
log(f"multiplier heights <= {h} digits")
log("STEP 11 PASS: Psi2(1/2) = (3/32)pi^2 - (9/8)Li2(1/3) + (3/4)log2log3 - (9/16)log^2 3 PROVED")
json.dump({"cols": [list(w) for w in cols], "prow": prow,
           "y": [str(yy) for yy in y],
           "statement": "Psi2(1/2) = (3/32)pi^2 - (9/8)Li2(1/3) + (3/4)log2 log3 - (9/16) log^2 3"},
          open("cert_psi2.json", "w"))
with open("progress.log", "a") as f:
    f.write("step_11: BONUS exact cert for Psi2 closed form (28 rows, residual 0 on 36 words) -> cert_psi2.json\n")
