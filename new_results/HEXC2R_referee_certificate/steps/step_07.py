"""step_07.py -- referee: rebuild T2/T3/T4 with MY OWN conversion AND MY OWN frame choices.
(1) exact-compare with hexc2_common targets (differences expected only from Li-atom frames);
(2) numeric: my targets vanish at 65 dps with MY evaluator;
(3) the differences my_T - their_T lie in the proved-relation span: EXACT (fmpq rref) at
    weights 2,3; mod 3 primes at weight 4;
(4) resolve the sign of H_{0,-1} in the printed ID2."""
import sys, json, re, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC2_level6_certificate")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
import mpmath as mp
from fractions import Fraction as Fr
import reflib as R
import hexc2_common as HC
import cmzv6
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
def my_fit_vec(name, weight):
    D = json.load(open(G1 + "/" + name))
    out = {}
    for key, val in D:
        m = re.match(r"c\^(-?\d+) H\[(.*)\]", key)
        if key == "1": e, w = 0, ()
        elif key == "c": e, w = 1, ()
        else:
            e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
        if len(w) != weight: continue
        v = {(): Fr(1)} if not w else R.my_hpl_to_cmzv(w)
        out = R.vadd((out, 1), (v, Fr(val) * Fr(1, 2) ** e))
    return out

def prod(*vecs):
    out = {(): Fr(1)}
    for v in vecs: out = R.mul(out, v)
    return out
def Li(n, z): return R.my_li_atom(n, Fr(z))   # MY frames

def my_F4_closed():
    L2, L3 = R.my_atom_log2(), R.my_atom_log3()
    Z3, PI2, PI4 = R.my_atom_zeta3(), {w: -c for w, c in R.my_atom_pi2().items()}, R.my_atom_pi4()
    Li21, Li31, Li41 = Li(2, '1/3'), Li(3, '1/3'), Li(4, '1/3')
    Li3m, Li4m = Li(3, '-1/2'), Li(4, '-1/2')
    out = {}
    def T(vec, coef):
        nonlocal out; out = R.vadd((out, 1), (vec, Fr(coef)))
    T(Li(4, -2), -21); T(Li41, Fr(-65, 8)); T(Li(4, '2/3'), Fr(13, 2)); T(Li(4, '1/2'), Fr(137, 6))
    T(Li41, Fr(39, 8)); T(Li(4, '-1/3'), Fr(39, 8))
    T(prod(Li31, L2), Fr(15, 2)); T(prod(Li3m, L2), Fr(37, 2)); T(prod(Li3m, L3), Fr(-5, 4))
    T(prod(Z3, L2), Fr(119, 16)); T(prod(Z3, L3), Fr(-37, 32))
    T(prod(Li21, Li21), 1); T(prod(Li21, L2, L2), Fr(-23, 4)); T(prod(Li21, L2, L3), Fr(7, 4))
    T(prod(Li21, L3, L3), 1); T(prod(PI2, Li21), Fr(-1, 6))
    T(prod(PI2, L2, L2), Fr(-455, 144)); T(prod(PI2, L2, L3), Fr(9, 8)); T(prod(PI2, L3, L3), Fr(-7, 32))
    T(PI4, Fr(-949, 2160))
    T(prod(L2, L2, L2, L2), Fr(-211, 144)); T(prod(L2, L2, L2, L3), Fr(101, 24))
    T(prod(L2, L2, L3, L3), Fr(-15, 8)); T(prod(L2, L3, L3, L3), Fr(-35, 24))
    T(prod(L3, L3, L3, L3), Fr(113, 192))
    return out

def my_W3_closed():
    L2, L3 = R.my_atom_log2(), R.my_atom_log3()
    Z3, PI2 = R.my_atom_zeta3(), {w: -c for w, c in R.my_atom_pi2().items()}
    Li21, Li31, Li3m = Li(2, '1/3'), Li(3, '1/3'), Li(3, '-1/2')
    out = {}
    def T(vec, coef):
        nonlocal out; out = R.vadd((out, 1), (vec, Fr(coef)))
    T(Z3, Fr(-15, 4)); T(Li31, -3); T(Li3m, -8)
    T(prod(PI2, L2), Fr(-1, 12)); T(prod(Li21, L2), 4); T(prod(Li21, L3), Fr(-3, 2))
    T(prod(L2, L2, L2), Fr(1, 3)); T(prod(L2, L2, L3), Fr(-5, 2)); T(prod(L2, L3, L3), 2)
    T(prod(L3, L3, L3), Fr(-1, 4))
    return out

def my_W2_closed():
    L2, L32 = R.my_atom_log2(), R.my_atom_log32()
    return R.vadd((prod(L2, L2), Fr(1, 2)), (Li(2, '1/3'), 1), (prod(L32, L32), Fr(1, 2)))

def my_target2():
    return R.vadd((my_fit_vec("hpl_fit_k3_(0,0,0).json", 2), 1),
                  (R.my_atom_zeta2(), Fr(1, 2)), (my_W2_closed(), -1))
def my_target3():
    R33 = my_fit_vec("hpl_fit_k3_(0,0,0).json", 3)
    Psi1 = my_fit_vec("hpl_fit_k3_(1,0,0).json", 1)
    lhs = R.vadd((R33, 1), (R.mul(R.my_atom_zeta2(), Psi1), 1), (R.my_atom_zeta3(), Fr(-1, 2)))
    return R.vadd((lhs, 1), (my_W3_closed(), -1))
def my_target4():
    D = json.load(open(G1 + "/hpl_fit_k3_(0,0,0).json"))
    hpl = {}
    for key, val in D:
        m = re.match(r"c\^(-?\d+) H\[(.*)\]", key)
        if not m: continue
        e = int(m.group(1))
        w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
        if len(w) != 4: continue
        hpl = R.vadd((hpl, 1), (R.my_hpl_to_cmzv(w), Fr(val) * Fr(1, 2) ** e))
    return R.vadd((hpl, 1), (my_F4_closed(), -1))

T_my = {2: my_target2(), 3: my_target3(), 4: my_target4()}
T_hc = {2: HC.target2(), 3: HC.target3(), 4: HC.target4()}
DIFF = {}
for n in (2, 3, 4):
    same = T_my[n] == T_hc[n]
    DIFF[n] = R.vadd((T_my[n], 1), (T_hc[n], -1))
    log(f"(1) my target{n} == theirs: {same} (supports mine {len(T_my[n])}, theirs {len(T_hc[n])}, diff {len(DIFF[n])})")

# (2) numeric vanishing of MY targets at 65 dps
mp.mp.dps = 65
for n in (2, 3, 4):
    v = R.cmzv_combo(T_my[n])
    log(f"(2) |my target{n} value| = {mp.nstr(abs(v), 3)}  (support {len(T_my[n])})")
    assert abs(v) < mp.mpf("1e-55")

# (3) differences lie in the proved span
from flint import fmpq_mat, fmpq
import numpy as np
lower = {}
def rank_modp(int_rows, ncols, p):
    """my own GF(p) Gaussian elimination (numpy int64): rank of the row set."""
    M = np.array(int_rows, dtype=np.int64) % p
    nr, rank, j = M.shape[0], 0, 0
    for j in range(ncols):
        piv = -1
        for i in range(rank, nr):
            if M[i, j]: piv = i; break
        if piv < 0: continue
        if piv != rank:
            M[[rank, piv]] = M[[piv, rank]]
        inv = pow(int(M[rank, j]), p - 2, p)
        M[rank] = (M[rank] * inv) % p
        # eliminate below only (row echelon, not reduced)
        for i in range(rank + 1, nr):
            if M[i, j]:
                M[i] = (M[i] - M[i, j] * M[rank]) % p
        rank += 1
        if rank == nr: break
    return rank

for n in (1, 2, 3, 4):
    rows, counts = cmzv6.make_relations(n, lower)
    lower[n] = rows
    if n < 2: continue
    cols = list(cmzv6.words(n))
    ids = {w: j for j, w in enumerate(cols)}
    D = DIFF[n]
    if not D:
        log(f"(3) weight {n}: targets identical, nothing to do")
        continue
    if n <= 3:
        # exact membership: rank(rows^T) vs rank(rows^T | D) over Q
        m = len(cols)
        A = fmpq_mat(m, len(rows))
        B = fmpq_mat(m, len(rows) + 1)
        for j, w in enumerate(cols):
            B[j, len(rows)] = fmpq(str(D.get(w, Fr(0))))
        for i, row in enumerate(rows):
            for w, c in row:
                if w in ids:
                    A[ids[w], i] = c
                    B[ids[w], i] = c
        r1 = A.rank(); r2 = B.rank()
        log(f"(3) weight {n}: rank(rows)={r1}, rank(rows+D)={r2}  "
            f"{'PASS: difference in span' if r1 == r2 else 'FAIL'}")
        assert r1 == r2
    else:
        int_rows = [[rowdict.get(w, 0) for w in cols] for rowdict in
                    ({w: c for w, c in row} for row in rows)]
        dvec = [D.get(w, Fr(0)) for w in cols]
        # clear denominators of D
        from math import lcm
        den = lcm(*[c.denominator for c in D.values()])
        dvec_i = [int(v * den) for v in dvec]
        for p in (1000000007, 1000000009, 998244353):
            r1 = rank_modp(int_rows, len(cols), p)
            r2 = rank_modp(int_rows + [dvec_i], len(cols), p)
            log(f"(3) weight {n} mod {p}: rank {r1} vs {r2}  "
                f"{'PASS: in span' if r1 == r2 else 'FAIL'}")
            assert r1 == r2

# (4) ID2 sign
mp.mp.dps = 50
H01 = R.hpl_value((0, 1), mp.mpf(1)/2)
H0m = R.hpl_value((0, -1), mp.mpf(1)/2)
rhs = mp.log(2)**2/2 + mp.polylog(2, mp.mpf(1)/3) + mp.log(mp.mpf(3)/2)**2/2
vA = -H01 + H0m + mp.zeta(2)/2
vB = -H01 - H0m + mp.zeta(2)/2
log(f"(4) fit-signed  -H01 + H0-1 + zeta(2)/2 - RHS = {mp.nstr(abs(vA-rhs),3)}")
log(f"    printed     -H01 - H0-1 + zeta(2)/2 - RHS = {mp.nstr(abs(vB-rhs),3)} (HEXC2 result.md ID2 line)")
assert abs(vA - rhs) < mp.mpf("1e-45") and abs(vB - rhs) > mp.mpf("1e-3")
log("    => certified ID2 has +H_{0,-1}; HEXC2 result.md display has a sign typo (harmless)")
log("STEP 07 PASS")
with open("progress.log", "a") as f:
    f.write("step_07: my-frame targets vanish numerically; diffs in proved span (exact w2/3, mod p w4); ID2 display typo\n")
