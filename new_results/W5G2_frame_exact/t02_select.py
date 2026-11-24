"""t02 (mod p): quotient map pi = (Fold_s(.) . N_s)_s from t01 (kernel = Fold(Et) mod p).
Candidates F' (all proved, literal): the 288 frame-w5 rows, the 32520 frame-234 lifts, and the u=7 images of the
frame-234 lifts.  Select a basis of pi(F') (greedy, in that order), solve pi(v) = sum a_j pi(F_j) mod p for the 74
targets v = u.f_i, and try rational reconstruction of the a_j from this one prime (a guess only; exactness comes later).
Writes sel.pkl."""
from g2common import *
import numpy as np, scipy.sparse as sps
from w5c2_fold import Ctx, sectorize_raw
from flint import nmod_mat
p = P1; LOG = 't02.log'
ctx = Ctx(); WW = ctx.W5; IDS = ctx.IDS
Q = {t: np.load(HERE + f'/quot_{t}.npz') for t in ('pp', 'pm', 'mp', 'mm')}
for t in Q: log(f'{t}: rank {int(Q[t]["rank"])} nullity {Q[t]["N"].shape[1]}', LOG)
def modq(c): c = F(c); return c.numerator % p * pow(c.denominator % p, -1, p) % p
def pi(vecs):
    blocks = []
    for t in ('pp', 'pm', 'mp', 'mm'):
        sec = SEC[t]; NC = len(ctx.sec_cols[sec])
        rr, cc, vv = [], [], []
        for i, v in enumerate(vecs):
            d = sectorize_raw([(w, modq(c)) for w, c in v.items()], sec, ctx)
            for k, a in d.items(): rr.append(i); cc.append(k); vv.append(a % p)
        M = sps.csr_matrix((np.array(vv, np.int64), (rr, cc)), shape=(len(vecs), NC))
        N = Q[t]['N']
        # exact mod-p product without overflow: split N into 12-bit halves
        lo = N & 4095; hi = N >> 12
        A = (M @ lo) % p; B = (M @ hi) % p
        blocks.append((A + B * 4096) % p)
    return np.hstack(blocks) % p
def act_id(d, u): return {IDS[act(WW[w], u)]: c for w, c in d.items()}
VW5 = load_framew5()
cand = []; labels = []
for i, row in enumerate(VW5['rows']):
    cand.append({IDS[tuple(w)]: c for w, c in row}); labels.append(('framew5', i, 1))
f234 = pickle.load(open(W5C1 + '/rels_frame234.pkl', 'rb'))
for i, (lab, r) in enumerate(f234):
    cand.append(dict(r)); labels.append(('frame234', i, 1))
for i, (lab, r) in enumerate(f234):
    cand.append(act_id(dict(r), 7)); labels.append(('frame234', i, 7))
t0 = time.time()
PI = pi(cand)
log(f'pi(F\'): {len(cand)} candidates, quotient dim {PI.shape[1]} ({time.time()-t0:.0f}s)', LOG)
# greedy basis: rref of the transpose, pivot columns = first independent candidates
nz = np.nonzero(PI.any(axis=1))[0]
Mt = nmod_mat(PI.shape[1], len(nz), [int(x) for x in PI[nz].T.ravel()], p)
rr, rk = Mt.rref()
piv = []; j = 0
for i in range(rk):
    while rr[i, j] == 0: j += 1
    piv.append(int(nz[j])); j += 1
log(f'rank pi(F\') = {rk}; basis families: ' + str(__import__("collections").Counter(labels[k][0] + str(labels[k][2]) for k in piv)), LOG)
Bas = PI[piv]                                   # rk x 408
tv = [{IDS[w]: c for w, c in target_vec(VW5, i, u).items()} for i, u in TARGETS]
PT = pi(tv)
# solve a.Bas = pt for each target: [Bas^T | PT^T]
Aug = nmod_mat(PI.shape[1], rk + len(TARGETS), [int(x) for x in np.hstack([Bas.T, PT.T]).ravel()], p)
ra, rka = Aug.rref()
log(f'rank [basis | targets] = {rka} (basis rank {rk}): all targets in span mod p: {rka == rk}', LOG)
assert rka == rk
# coefficients: rref of [Bas^T | PT^T]; since Bas^T has full column rank rk, first rk pivots are columns 0..rk-1
a = np.array([[int(ra[i, rk + t]) for t in range(len(TARGETS))] for i in range(rk)], np.int64)   # rk x 74
def ratrec(x, m=p):
    import math
    N = math.isqrt(m // 2); r0, r1, t0_, t1 = m, x % m, 0, 1
    while r1 > N:
        q = r0 // r1; r0, r1 = r1, r0 - q * r1; t0_, t1 = t1, t0_ - q * t1
    if t1 == 0 or abs(t1) > N or math.gcd(r1, abs(t1)) != 1: return None
    return F(r1, t1)
rec = [[ratrec(int(a[i, t])) for t in range(len(TARGETS))] for i in range(rk)]
nrec = sum(1 for i in range(rk) for t in range(len(TARGETS)) if rec[i][t] is not None)
used = sorted({i for i in range(rk) for t in range(len(TARGETS)) if a[i, t]})
log(f'a mod p: nonzero basis elements used by some target: {len(used)}/{rk}; entries reconstructible with |num|,|den|<=2896: {nrec}/{rk*len(TARGETS)}', LOG)
for t, (i, u) in enumerate(TARGETS):
    nzr = [(labels[piv[k]], rec[k][t]) for k in range(rk) if a[k, t]]
    log(f'  u{u}.f{i}: {len(nzr)} frame terms; ' + ', '.join(f'{l}:{r}' for l, r in nzr[:12]) + (' ...' if len(nzr) > 12 else ''), LOG)
pickle.dump(dict(piv=[labels[k] for k in piv], a_modp=a, rec=rec, rank=rk), open(HERE + '/sel.pkl', 'wb'))
log('saved sel.pkl', LOG)
