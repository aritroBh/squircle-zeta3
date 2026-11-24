import os, pickle, numpy as np, flint, sys
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
sys.path.insert(0,'.')
d = pickle.load(open('work/step01_pp.pkl','rb'))
S, C = d['S'], d['C']
feed = pickle.load(open('work/feed_pp.pkl','rb'))
p = 16777213
Cset = sorted(C); cpos = {c:k for k,c in enumerate(Cset)}
r = len(S)
A = np.zeros((r, r), np.int64)
for j, si in enumerate(S):
    for c, a in feed[si][2]:
        k = cpos.get(c)
        if k is not None:
            A[k, j] = a
print('A built, nnz', (A!=0).sum(), 'maxA', abs(A).max())
# rank via flint nmod
M = flint.nmod_mat((A % p).astype(object).tolist(), p)
print('flint nmod rank:', M.rank(), 'of', r)