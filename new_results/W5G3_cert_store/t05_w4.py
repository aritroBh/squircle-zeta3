"""t05 (EXACT): every weight-4 vector W (from t04: v = proved + e sh W) is an exact rational combination of proved
weight-4 relations  G = A3h5's 5042 base rows T0[S] (std 4810, lift1 223, lift2old 7, lift3old 2; Theorem W4R's
inputs)  plus the weight-4 lifts of A3h4's validated frame rows of weights 2, 3, 4 (proved like the 28 frame lifts
of Theorem W4R; each one used is regenerated in domain by t06).  Then e sh W = sum t_k (e sh G_k) and every
e sh G_k is a relation (shuffle of convergent iterated integrals: value (i pi/2) * 0).
Method: independent generators + pivot columns mod p, Dixon on the square system, exact check on all 5184 words."""
from g2lin import *
import certio
import numpy as np, scipy.sparse as sps, collections, math
import dixon2
from flint import nmod_mat
LOG = 't05.log'; p = P1
t0 = time.time()
WV = pickle.load(open(HERE + '/wvecs.pkl', 'rb'))
keys = [k for k in TARGETS if WV.get(k)]
log(f'W vectors: {len(WV)}, nonzero: {len(keys)}', LOG)
B = bases(); W4 = words(4); ids = {w: i for i, w in enumerate(W4)}
gens = list(B[4]); labels = [('A3h5_S', i) for i in range(len(B[4]))]
for n in (2, 3, 4):
    for i, row in enumerate(pickle.load(open(A3H4 + f'/validated_w{n}.pickle', 'rb'))['rows']):
        for w in words(4 - n):
            gens.append({k: F(v) for k, v in mul(dict(row), {w: 1}).items()}); labels.append(('A3h4_frame', n, i, w))
gl = [math.lcm(*[c.denominator for c in g.values()]) if g else 1 for g in gens]
G = [{ids[w]: int(c * gl[a]) for w, c in g.items() if c} for a, g in enumerate(gens)]
log(f'generators {len(G)} ({time.time()-t0:.0f}s)', LOG)
Mm = nmod_mat(len(G), len(W4), p)
for a, g in enumerate(G):
    for c, v in g.items(): Mm[a, c] = v % p
tr, rk = Mm.transpose().rref()
I = []; j = 0
for i in range(rk):
    while tr[i, j] == 0: j += 1
    I.append(j); j += 1
MI = nmod_mat(rk, len(W4), p)
for a, i in enumerate(I):
    for c, v in G[i].items(): MI[a, c] = v % p
rr, rk2 = MI.rref(); assert rk2 == rk
C = []; j = 0
for i in range(rk):
    while rr[i, j] == 0: j += 1
    C.append(j); j += 1
cpos = {c: k for k, c in enumerate(C)}
log(f'rank {rk}; independent generators chosen; frame lifts among them: {sum(1 for i in I if labels[i][0] != "A3h5_S")} ({time.time()-t0:.0f}s)', LOG)
rows_, cols_, vals_ = [], [], []
for a, i in enumerate(I):
    for c, v in G[i].items():
        if c in cpos: rows_.append(cpos[c]); cols_.append(a); vals_.append(v)
M = sps.csr_matrix((np.array(vals_, np.int64), (rows_, cols_)), shape=(rk, rk))
# membership mod p first (rank test), then exact
Ls = []; Bm = np.zeros((rk, len(keys)), dtype=object); Bm[:] = 0
for t, key in enumerate(keys):
    Wv = WV[key]; L = 1
    for c in Wv.values(): L = math.lcm(L, F(c).denominator)
    Ls.append(L)
    for w, c in Wv.items():
        if ids[w] in cpos: Bm[cpos[ids[w]], t] = int(F(c) * L)
Af = dixon2.modinv_dense(M, p, None, HERE + '/' + LOG)
n, d, steps = dixon2.dixon_solve(M, Bm, p, Af, logf=HERE + '/' + LOG)
summary = {}; used_frames = set(); allok = True
for t, key in enumerate(keys):
    coef = {I[a]: F(int(n[a, t]), d[t] * Ls[t]) * gl[I[a]] for a in range(rk) if n[a, t]}
    acc = collections.defaultdict(lambda: F(0))
    for i, c in coef.items():
        for w, v in gens[i].items(): acc[w] += c * v
    Wv = WV[key]
    bad = sum(1 for w in set(acc) | set(Wv) if acc.get(w, F(0)) != F(Wv.get(w, 0)))
    fr = sorted({labels[i][1:3] for i in coef if labels[i][0] == 'A3h4_frame'})
    used_frames |= set(fr)
    hd = max(len(str(abs(c.numerator))) + len(str(c.denominator)) for c in coef.values())
    summary[key] = dict(n_gens=len(coef), n_frame_lifts=sum(1 for i in coef if labels[i][0] == 'A3h4_frame'),
                        frame_rows=fr, max_digits=hd, mismatches=bad)
    allok &= bad == 0
    # ---- W5G3: append section t (coefficients t_k with the identity of each generator G_k) ----
    tl = []
    for gi in sorted(coef):
        lb = labels[gi]
        if lb[0] == 'A3h5_S': tl.append(['A3h5_S', lb[1], certio.fstr(coef[gi])])
        else: tl.append(['A3h4_frame', lb[1], lb[2], list(lb[3]), certio.fstr(coef[gi])])
    certio.write_lines(certio.path(key[0], key[1], 'core'), iter([certio.jl(dict(section='t', n_gens=len(tl), mismatches_w4=bad, terms=tl))]), 'ab')
    log(f'W(u{key[1]}.f{key[0]}): |supp| {len(Wv)}; = exact combination of {len(coef)} proved weight-4 rows '
        f'({summary[key]["n_frame_lifts"]} frame lifts, frame rows {fr}); EXACT check on 5184 words: mismatches {bad}', LOG)
for key in TARGETS:
    if key not in summary: summary[key] = dict(n_gens=0, n_frame_lifts=0, frame_rows=[], max_digits=0, mismatches=0)
pickle.dump(dict(summary=summary, frame_rows_used=sorted(used_frames)), open(HERE + '/w4cert_summary.pkl', 'wb'))
log(f'ALL W in the proved weight-4 span (exact): {allok} ({time.time()-t0:.0f}s)', LOG)
assert allok
