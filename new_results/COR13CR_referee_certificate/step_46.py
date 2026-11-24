import time, sys, json
sys.dont_write_bytecode = True
base = '<project>/research/08_theory'
sys.path.insert(0, base + '/COR13C_third_corner_certificate')
sys.path.insert(0, base + '/HEXC1_hexagon_certificate')
import mpmath as mp
from fractions import Fraction as Fr
import reflib13 as R13
import reflib as R
import cor13c_common as HC
import cmzv6
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

# (1) with HC's frames, my target4 == HC target4 exactly: rebuild target4 with HC frames
def li_atom_hcframe(n, z):
    frame = HC.cross_table()[Fr(z)][0]
    frame = tuple(None if p is HC.INF else p for p in frame)
    return R13.my13_li_atom(n, Fr(z), frame=frame)
def atom_vec_hcf(name):
    if name == "zeta3": return R13.a_zeta3()
    if name == "pi4": return R13.a_pi4()
    m = re.fullmatch(r"Li(\d)\((-?\d/\d|-?\d+)\)", name)
    if m:
        n, zz = int(m.group(1)), Fr(m.group(2))
        if zz == Fr(1,9):
            return R13.smul(Fr(8), R.vadd((li_atom_hcframe(4,"1/3"),1),(li_atom_hcframe(4,"-1/3"),1)))
        return li_atom_hcframe(n, zz)
    facs = []
    for f in name.split("*"):
        mm = re.fullmatch(r"(Li\d\(-?\d/\d\)|zeta3|pi2|log2|log3)(?:\^(\d+))?", f)
        base_, pw = mm.group(1), int(mm.group(2) or 1)
        if base_ == "pi2": v = R13.smul(Fr(-1), R13.a_pi2())
        elif base_ == "log2": v = R13.a_log2()
        elif base_ == "log3": v = R13.a_log3()
        elif base_ == "zeta3": v = R13.a_zeta3()
        else:
            n2, zz = int(base_[2]), Fr(base_[4:-1]); v = li_atom_hcframe(n2, zz)
        facs += [v]*pw
    return R13.prod(*facs)
import re
W4_hcf = {}
Rf = json.load(open(base + "/COR7_classical_corners/c13_F4_relation.json"))
for cf, bn in zip(Rf["coeffs"], Rf["basis"]):
    W4_hcf = R.vadd((W4_hcf, 1), (atom_vec_hcf(bn), Fr(-cf, Rf["a0"])))
T4_hcf = R.vadd((R13.R3_13(4), 1), (R.mul(R13.a_zeta2(), R13.Psi_13(2)), 1),
                (R.mul(R13.a_zeta3(), R13.a_log2()), Fr(2,3)), (W4_hcf, -1))
T4_hc = {w: Fr(c) for w, c in HC.target4().items()}
log(f"(1) my target4 with HC frames == cor13c target4 exactly: {T4_hcf == T4_hc} (supports {len(T4_hcf)}/{len(T4_hc)})")
assert T4_hcf == T4_hc

# (2) MY-frame targets: numeric vanishing at 65 dps with my evaluator
mp.mp.dps = 65
T_my = {2: R13.my_target2(), 3: R13.my_target3(), 4: R13.my_target4()}
for n in (2, 3, 4):
    v = R.cmzv_combo(T_my[n])
    log(f"(2) |my-frame target{n} value| = {mp.nstr(abs(v),3)} (support {len(T_my[n])})")
    assert abs(v) < mp.mpf('1e-55')

# (3) differences (my frames - their frames) lie in the proved span
T_hc = {2: {w: Fr(c) for w, c in HC.target2().items()},
        3: {w: Fr(c) for w, c in HC.target3().items()},
        4: T4_hc}
DIFF = {n: R.vadd((T_my[n], 1), (T_hc[n], -1)) for n in (2,3,4)}
log(f"(3) difference supports: {[len(DIFF[n]) for n in (2,3,4)]}")
lower = {}
import numpy as np
def rank_modp(int_rows, ncols, p):
    M = np.array(int_rows, dtype=np.int64) % p
    nr, rank = M.shape[0], 0
    for j in range(ncols):
        if rank >= nr: break
        piv = -1
        for i in range(rank, nr):
            if M[i, j]: piv = i; break
        if piv < 0: continue
        if piv != rank: M[[rank, piv]] = M[[piv, rank]]
        inv = pow(int(M[rank, j]), p-2, p)
        M[rank] = (M[rank]*inv) % p
        for i in range(rank+1, nr):
            if M[i, j]: M[i] = (M[i] - M[i, j]*M[rank]) % p
        rank += 1
    return rank
from flint import fmpq_mat
for n in (1, 2, 3, 4):
    rows, counts = cmzv6.make_relations(n, lower)
    lower[n] = rows
    if n < 2: continue
    cols = list(cmzv6.words(n))
    ids = {w: j for j, w in enumerate(cols)}
    Dv = DIFF[n]
    if not Dv:
        log(f"(3) weight {n}: targets identical, nothing to check"); continue
    if n <= 3:
        m = len(cols)
        A = fmpq_mat(m, len(rows)); B = fmpq_mat(m, len(rows)+1)
        for j, w in enumerate(cols):
            B[j, len(rows)] = fmpq_mat(1,1,[str(Dv.get(w, Fr(0)))])[0,0]
        for i, row in enumerate(rows):
            for w, c in row:
                if w in ids:
                    A[ids[w], i] = c; B[ids[w], i] = c
        r1, r2 = A.rank(), B.rank()
        log(f"(3) weight {n}: rank(rows)={r1} rank(rows+D)={r2} {'PASS in span' if r1==r2 else 'FAIL'}")
        assert r1 == r2
    else:
        int_rows = [[dict(row).get(w, 0) for w in cols] for row in rows]
        from math import lcm
        den = lcm(*[c.denominator for c in Dv.values()])
        dvec = [int(Dv.get(w, Fr(0))*den) for w in cols]
        for p in (1000000007, 1000000009, 998244353):
            r1 = rank_modp(int_rows, len(cols), p)
            r2 = rank_modp(int_rows + [dvec], len(cols), p)
            log(f"(3) weight {n} mod {p}: rank {r1} vs {r2} {'PASS in span' if r1==r2 else 'FAIL'}")
            assert r1 == r2
log("STEP PASS: frame independence (my frames vanish numerically; differences in proved span)")
with open('progress.log','a') as f:
    f.write("referee step(frames): my-frame targets vanish at 65dps; frame differences in proved span (exact w2/3, mod 3 primes w4); HC-frame targets == theirs exactly\n")
