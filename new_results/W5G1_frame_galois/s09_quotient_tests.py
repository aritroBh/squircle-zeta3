"""s09 (route C2, mod p, uses s08's quotient maps pi = (Fold_s(.) . N_s)_s, kernel = Et).
P' = the 288 frame-w5 rows (u=7-closed) + all frame-234 lifts (rels_frame234.pkl) + their u=7 images.
Tests: rank pi(P'); for every frame (row,u) pair used by W5C2 (u=3,5): is pi(u.f) in span pi(P')?;
is pi(X_rat), pi(X_sq) in span pi(P')?  (X from s05: D = justified + e sh Z + X.)"""
from common import *
sys.set_int_max_str_digits(0)
import numpy as np, scipy.sparse as sps
sys.path.insert(0, W5C2)
from w5c2_fold import Ctx, sectorize_raw
from w5c2_lib import P1
from flint import nmod_mat
p = P1
ctx = Ctx(); W5 = ctx.W5; IDS = ctx.IDS
Q = {t: np.load(HERE + f'/quot_{t}.npz') for t in ('pp', 'pm', 'mp', 'mm')}
for t in Q: log(f'{t}: rank {int(Q[t]["rank"])} nullity {Q[t]["N"].shape[1]}', 's09.log')
def modq(c): c = F(c); return c.numerator % p * pow(c.denominator % p, -1, p) % p
def pi(vecs):
    """vecs: list of dict word_id -> coeff. returns (len, sum q) int64 mod p."""
    blocks = []
    for t in ('pp', 'pm', 'mp', 'mm'):
        sec = SEC[t]; NC = len(ctx.sec_cols[sec])
        rr, cc, vv = [], [], []
        for i, v in enumerate(vecs):
            d = sectorize_raw([(w, modq(c)) for w, c in v.items()], sec, ctx)
            for k, a in d.items(): rr.append(i); cc.append(k); vv.append(a % p)
        M = sps.csr_matrix((np.array(vv, np.int64), (rr, cc)), shape=(len(vecs), NC))
        blocks.append((M @ Q[t]['N']) % p)
    return np.hstack(blocks) % p
def rank(A):
    if A.shape[0] == 0: return 0
    return nmod_mat(A.shape[0], A.shape[1], [int(x) for x in A.ravel()], p).rank()
VW5 = pickle.load(open(W5C1 + '/validated_w5.pkl', 'rb'))
Fr = [{IDS[tuple(w)]: c for w, c in row} for row in VW5['rows']]
f234 = [dict(r) for lab, r in pickle.load(open(W5C1 + '/rels_frame234.pkl', 'rb'))]
def act_id(d, u): return {IDS[act(W5[w], u)]: c for w, c in d.items()}
Pp = Fr + f234 + [act_id(d, 7) for d in f234]
t0 = time.time()
PI = pi(Pp)
# reduce P' image to a basis
B = nmod_mat(PI.shape[0], PI.shape[1], [int(x) for x in PI.ravel()], p)
rr, rk = B.rref()
basis = np.array([[int(rr[i, j]) for j in range(PI.shape[1])] for i in range(rk)], np.int64)
log(f'pi(P\'): {len(Pp)} vectors in quotient of dim {PI.shape[1]}; rank {rk}; frame-w5 alone rank {rank(PI[:len(Fr)])}; f234 part rank {rank(PI[len(Fr):])} ({time.time()-t0:.0f}s)', 's09.log')
sp = pickle.load(open(HERE + '/split.pkl', 'rb'))
pairs = sorted(sp['frame_pairs'])
tests = [act_id(Fr[i], u) for i, u in pairs]
names = [f'u{u}.f{i}' for i, u in pairs]
for c in ('rat', 'sq'):
    tests.append({IDS[w]: v for w, v in sp['X'][c].items() if v}); names.append('X_' + c)
PT = pi(tests)
res = {}
for nm, v in zip(names, PT):
    res[nm] = (rank(np.vstack([basis, v[None, :]])) == rk, bool(np.any(v)))
log('membership pi(v) in span pi(P\') [in_span, v_not_in_Et]: ' + str(res), 's09.log')
allF = np.vstack([PI[:len(Fr)], PT[:len(pairs)]])
log(f'rank pi(F)={rank(PI[:len(Fr)])}, rank pi(F + used u.f)={rank(allF)}, rank pi(P\' + used u.f)={rank(np.vstack([basis, PT[:len(pairs)]]))}', 's09.log')
pickle.dump(dict(res=res, rk=rk), open(HERE + '/s09.pkl', 'wb'))
