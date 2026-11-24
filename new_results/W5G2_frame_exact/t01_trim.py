"""t01 (mod p, one sector): fold the G-stable generator set
     Et = std (weight-5 double shuffle + distribution), lift1, I = e sh words(4), lift2old, lift3old, lift4old
into sector s, echelon mod P1 with pivot-row tracking (LeanEchelon, from W5C2), then the right null basis N_s
(back substitution) for the quotient map pi_s(b) = Fold_s(b).N_s  (kernel = Fold_s(Et) mod p).
Feed order: std, lift1, I (shortest first), then lift2old/3old/4old (shortest first) so that S uses as few
higher lifts as possible (each needs an exact decomposition of its base row's image later).
Writes trim_<s>.pkl (S rows with provenance, C, free cols) and quot_<s>.npz.   Usage: t01_trim.py pp|pm|mp|mm"""
from g2common import *
import numpy as np
from w5c2_fold import Ctx, sectorize_raw, norm
from w5c2_lib import LeanEchelon
tag = sys.argv[1]; sec = SEC[tag]; p = P1; LOG = f't01_{tag}.log'
t0 = time.time()
ctx = Ctx(); ctx.check_vs_meta()
W = ctx.W5; IDS = ctx.IDS; NC = len(ctx.sec_cols[sec])
gens = []   # (fam, idx, unfolded row as tuple((word_id,int)))
for f in ('std', 'lift1'):
    gens += [(f, i, r) for i, (lab, r) in enumerate(pickle.load(open(W5C1 + f'/rels_{f}.pkl', 'rb')))]
for i, w in enumerate(words(4)):
    d = mul(E1, {w: 1})
    gens.append(('I', i, tuple(sorted((IDS[x], int(c)) for x, c in d.items() if c))))
for f in ('lift2old', 'lift3old', 'lift4old'):
    gens += [(f, i, r) for i, (lab, r) in enumerate(pickle.load(open(W5C1 + f'/rels_{f}.pkl', 'rb')))]
log(f'{tag}: unfolded generators {len(gens)}, sector cols {NC} ({time.time()-t0:.0f}s)', LOG)
seen = set(); rows = []
for f, i, r in gens:
    d = sectorize_raw(r, sec, ctx)
    if not d: continue
    k, lam = norm(d)
    if k in seen: continue
    seen.add(k)
    grp = 0 if f in ('std', 'lift1', 'I') else 1
    rows.append((grp, len(k), f, i, lam, k))
rows.sort(key=lambda t: (t[0], t[1]))
log(f'{tag}: folded unique rows {len(rows)} ({time.time()-t0:.0f}s)', LOG)
E = LeanEchelon(NC, p)
for ri, (grp, ln, f, i, lam, k) in enumerate(rows):
    cols = np.array([c for c, _ in k], np.int64); vals = np.array([v % p for _, v in k], np.int64)
    buf = E.reduce(cols, vals)
    if buf.any(): E.add_pivot(buf, src=ri)
    if (ri + 1) % 10000 == 0:
        log(f'{tag}: row {ri+1}/{len(rows)} rank {E.r} ({time.time()-t0:.0f}s)', LOG)
r = E.r
pc = [int(c) for c in E.pivcol[:r]]
free = sorted(set(range(NC)) - set(pc))
q = len(free)
import collections
famc = collections.Counter(rows[s][2] for s in E.pivrow)
log(f'{tag}: rank {r}, nullity {q}; S families {dict(famc)} ({time.time()-t0:.0f}s)', LOG)
X = np.zeros((NC, q), np.int64)
for j, c in enumerate(free): X[c, j] = 1
for t in range(r - 1, -1, -1):
    idx = E.piv_idx[t].astype(np.int64); v = E.piv_val[t]; c = pc[t]
    m = idx != c
    X[c, :] = (-(v[m][:, None] * X[idx[m], :]).sum(axis=0)) % p
# verify null basis on ALL folded rows (mod p)
bad = 0
for ri in range(len(rows)):
    k = rows[ri][5]
    cc = np.array([c for c, _ in k]); vv = np.array([v % p for _, v in k], np.int64)
    s = (vv[:, None] * X[cc, :]).sum(axis=0) % p
    if s.any(): bad += 1
log(f'{tag}: null basis checked on all {len(rows)} folded rows mod p: bad={bad} ({time.time()-t0:.0f}s)', LOG)
assert bad == 0
np.savez_compressed(HERE + f'/quot_{tag}.npz', free=np.array(free), N=X, rank=r)
Srows = [rows[s][2:] for s in E.pivrow]          # (fam, idx, lam, key)
pickle.dump(dict(S=Srows, C=pc, free=free, rank=r, ncols=NC, p=p), open(HERE + f'/trim_{tag}.pkl', 'wb'), protocol=4)
log(f'{tag}: saved trim_{tag}.pkl, quot_{tag}.npz ({time.time()-t0:.0f}s)', LOG)
