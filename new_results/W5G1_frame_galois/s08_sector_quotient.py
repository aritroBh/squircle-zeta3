"""s08 (route C2, mod p): Galois-sector quotient of the weight-5 word space by the G-STABLE space
   Et = span( std rows, all lift rows (lift1, lift2old, lift3old, lift4old), I ),   I = e sh Q^{words(4)}.
Et is G-stable: std is literally closed; u.(b sh w) = (rho sh u.w) + e sh (y sh u.w) (s03/s04: every u.b lies in
P_n + e sh words(n-1)); u.(e sh w) = +-e sh u.w.  NOTE: I is not made of relations (e.V = i*pi/2); it is only used
to get a G-stable space; whatever I-part a membership needs is accounted for afterwards (route via Z in P_4).
For each sector s: fold all generators (W5C2's fold, exactly), sparse echelon mod p (W5C2 LeanEchelon), then a
basis N_s of the right null space (back substitution).  v in Et  <=>  Fold_s(v).N_s = 0 for all s.
Writes quot_<s>.npz (free cols, N_s) and logs ranks.  Usage: s08_sector_quotient.py [pp pm mp mm]"""
from common import *
import numpy as np
sys.path.insert(0, W5C2)
from w5c2_fold import Ctx, sectorize_raw, norm, SECTORS
from w5c2_lib import LeanEchelon, P1
p = P1
ctx = Ctx(); ctx.check_vs_meta()
W5 = ctx.W5; IDS = ctx.IDS
tags = sys.argv[1:] or ['pp', 'pm', 'mp', 'mm']
t0 = time.time()
unf = []
for f in ('std', 'lift1', 'lift2old', 'lift3old', 'lift4old'):
    unf += [r for lab, r in pickle.load(open(W5C1 + f'/rels_{f}.pkl', 'rb'))]
nE = len(unf)
for w in words(4):
    d = mul(E1, {w: 1})
    unf.append(tuple(sorted((IDS[x], int(c)) for x, c in d.items())))
log(f's08: unfolded generators {len(unf)} (engine {nE}, I {len(unf)-nE}) {time.time()-t0:.0f}s', 's08.log')
for tag in tags:
    sec = SEC[tag]; NC = len(ctx.sec_cols[sec]); ts = time.time()
    seen = set(); rows = []
    for r in unf:
        d = sectorize_raw(r, sec, ctx)
        if not d: continue
        k, lam = norm(d)
        if k in seen: continue
        seen.add(k); rows.append(k)
    rows.sort(key=len)
    log(f'{tag}: folded unique rows {len(rows)} cols {NC} ({time.time()-ts:.0f}s)', 's08.log')
    E = LeanEchelon(NC, p)
    for ri, k in enumerate(rows):
        cols = np.array([c for c, _ in k], np.int64); vals = np.array([v % p for _, v in k], np.int64)
        buf = E.reduce(cols, vals)
        if buf.any(): E.add_pivot(buf, src=ri)
        if (ri + 1) % 10000 == 0:
            log(f'{tag}: row {ri+1}/{len(rows)} rank {E.r} ({time.time()-ts:.0f}s)', 's08.log')
    r = E.r
    pc = [int(c) for c in E.pivcol[:r]]
    free = sorted(set(range(NC)) - set(pc))
    q = len(free)
    log(f'{tag}: rank {r}, nullity {q} ({time.time()-ts:.0f}s)', 's08.log')
    X = np.zeros((NC, q), np.int64)
    for j, c in enumerate(free): X[c, j] = 1
    for t in range(r - 1, -1, -1):
        idx = E.piv_idx[t].astype(np.int64); v = E.piv_val[t]; c = pc[t]
        m = idx != c
        # pivot value normalised to 1
        X[c, :] = (-(v[m][:, None] * X[idx[m], :]).sum(axis=0)) % p
    # verify on a sample of original rows
    rng = np.random.default_rng(1)
    for ri in rng.choice(len(rows), 300, replace=False):
        k = rows[ri]
        s = sum(int(v) * X[c, :] for c, v in k) % p
        assert not np.any(s), (tag, ri)
    np.savez_compressed(HERE + f'/quot_{tag}.npz', free=np.array(free), N=X, rank=r)
    log(f'{tag}: null basis verified on 300 sample rows; saved quot_{tag}.npz ({time.time()-ts:.0f}s)', 's08.log')
