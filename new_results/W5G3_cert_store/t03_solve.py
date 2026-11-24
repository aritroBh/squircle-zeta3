"""t03 (EXACT, one sector s): for every RHS b in
     {f_i : the 37 open frame-w5 rows}  u  {the frame basis F_j chosen in t02}
solve y_b . R_{S,C} = Fold_s(b)_C exactly over Q (Dixon, dixon2.py), and compute the exact normal form
     NF_s(b) = Fold_s(b)_free - y_b . R_{S,free}        (zero  <=>  Fold_s(b) in Fold_s(span S)).
Note Fold_s(u.f) = chi_s(u) Fold_s(f), so the 74 targets need only the 37 rows f_i as right-hand sides.
R = normalised folded generator keys of t01 (S rows).  Writes ysol_<s>.pkl.   Usage: t03_solve.py pp|pm|mp|mm"""
from g2common import *
import numpy as np, scipy.sparse as sps
from w5c2_fold import Ctx, sectorize_raw
import dixon2
tag = sys.argv[1]; sec = SEC[tag]; p = P1; LOG = f't03_{tag}.log'
t0 = time.time()
ctx = Ctx(); WW = ctx.W5; IDS = ctx.IDS
tr = pickle.load(open(SRC + f'/trim_{tag}.pkl', 'rb'))
S, C, free = tr['S'], tr['C'], tr['free']
r = len(S); cpos = {c: k for k, c in enumerate(C)}; fpos = {c: k for k, c in enumerate(free)}
sel = pickle.load(open(SRC + '/sel.pkl', 'rb'))
VW5 = load_framew5()
f234 = None
def unfolded(lab):
    global f234
    fam, i, u = lab
    if fam == 'framew5': d = {IDS[tuple(w)]: c for w, c in VW5['rows'][i]}
    else:
        if f234 is None: f234 = pickle.load(open(W5C1 + '/rels_frame234.pkl', 'rb'))
        d = dict(f234[i][1])
    if u != 1: d = {IDS[act(WW[w], u)]: c for w, c in d.items()}
    return d
rhs_labels = [('framew5', i, 1) for i in OPEN_ROWS] + [l for l in sel['piv'] if l not in [('framew5', i, 1) for i in OPEN_ROWS]]
m = len(rhs_labels)
B = np.zeros((r, m), np.int64); BF = np.zeros((len(free), m), np.int64)
for j, lab in enumerate(rhs_labels):
    d = sectorize_raw(list(unfolded(lab).items()), sec, ctx)
    for k, a in d.items():
        if k in cpos: B[cpos[k], j] = a
        else: BF[fpos[k], j] = a
log(f'{tag}: r={r} free={len(free)} rhs={m} (37 open rows + {m-37} basis); maxB={int(abs(B).max())} ({time.time()-t0:.0f}s)', LOG)
rr, cc, vv = [], [], []; RF = {}
for i, (fam, idx, lam, key) in enumerate(S):
    for c, a in key:
        if c in cpos: rr.append(cpos[c]); cc.append(i); vv.append(a)
        else: RF.setdefault(i, []).append((fpos[c], a))
M = sps.csr_matrix((np.array(vv, np.int64), (rr, cc)), shape=(r, r))      # M = R_{S,C}^T
Af = dixon2.modinv_dense(M, p, HERE + f'/ainv_{tag}.npy', HERE + '/' + LOG)
n, d, steps = dixon2.dixon_solve(M, B, p, Af, logf=HERE + '/' + LOG)
del Af
# exact normal forms
NF = []
for j in range(m):
    acc = [0] * len(free)
    for i, lst in RF.items():
        x = n[i, j]
        if x:
            for k, a in lst: acc[k] += a * x
    NF.append({free[k]: F(int(BF[k, j]) * d[j] - acc[k], d[j]) for k in range(len(free)) if int(BF[k, j]) * d[j] != acc[k]})
hn = max(len(str(abs(int(x)))) for x in n.ravel()); hd = max(len(str(x)) for x in d)
log(f'{tag}: exact y for {m} rhs, Dixon steps {steps}; max numerator digits {hn}, max common-denominator digits {hd}; '
    f'NF nonzero for {sum(1 for x in NF if x)}/{m} rhs ({time.time()-t0:.0f}s)', LOG)
pickle.dump(dict(rhs=rhs_labels, n=n, d=d, NF=NF, S=S, C=C, free=free), open(HERE + f'/ysol_{tag}.pkl', 'wb'), protocol=4)
log(f'{tag}: saved ysol_{tag}.pkl ({time.time()-t0:.0f}s)', LOG)
