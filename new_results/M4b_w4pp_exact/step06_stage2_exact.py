"""step06: Stage 2 exact. G_X = W4pp_X - sum_k c_k^X B_k - sum_j z_j^X F_j  (exact Fractions);
solve y^T S = G against the 5042 independent T0 rows S (A3h5 stage2_rows_S.pkl) by exact fmpz solve
on pivot columns pp (rref_T0_1000000009.pkl); then EXACT Fraction verification on ALL 5184 columns.
Saves stage2_G.pkl, stage2_y.pkl."""
import sys, json, pickle, time, bisect
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/M4b_w4pp_exact')
import m4b_common as mc
from fractions import Fraction as F
from math import lcm
import flint

t0 = time.time()
W = pickle.load(open(str(mc.HERE/'stage1_w4pp_vectors.pkl'), 'rb'))
AT = pickle.load(open(str(mc.HERE/'stage2_atoms.pkl'), 'rb'))['atoms']
d = json.load(open(str(mc.HERE/'stage1_exact_cz.json')))
cr = [F(x) for x in d['c_rat']]; cs = [F(x) for x in d['c_sq']]
zr = [F(x) for x in d['z_rat']]; zs = [F(x) for x in d['z_sq']]
chosen = d['chosen_frames']
frames = mc.load_frame_rows()
for k, idx in enumerate(chosen):
    assert str(frames[idx][0]) == d['frame_labels'][k]
Bvecs = [v for _, _, v in AT]
Fvecs = [frames[i][1] for i in chosen]

def build(cc, zz):
    G = dict(W['rat' if cc is cr else 'sq'])
    for k, a in enumerate(cc):
        if not a: continue
        for w, c in Bvecs[k].items():
            G[w] = G.get(w, F(0)) - a*c
            if not G[w]: del G[w]
    for j, a in enumerate(zz):
        if not a: continue
        for w, c in Fvecs[j]:
            G[w] = G.get(w, F(0)) - a*c
            if not G[w]: del G[w]
    return G
Grat = build(cr, zr); Gsq = build(cs, zs)
cols = mc.words4(); cset = set(cols)
assert all(w in cset for w in Grat) and all(w in cset for w in Gsq)
Lrat = 1; Lsq = 1
for c in Grat.values(): Lrat = lcm(Lrat, c.denominator)
for c in Gsq.values(): Lsq = lcm(Lsq, c.denominator)
mc.log(f'G supports rat={len(Grat)} sq={len(Gsq)}; Lrat={Lrat} Lsq={Lsq} s={time.time()-t0:.1f}')
with open(str(mc.HERE/'stage2_G.pkl'), 'wb') as f:
    pickle.dump({'Grat': Grat, 'Gsq': Gsq, 'Lrat': Lrat, 'Lsq': Lsq}, f)

# system: rows S, pivot columns pp
rows = mc.load_T0()
S = pickle.load(open(str(mc.THEORY/'A3h5_exact_certificate/stage2_rows_S.pkl'), 'rb'))['S']
R = mc.load(mc.rref_path(mc.P1))
pp = R['pp']
assert len(S) == len(pp) == 5042
ids = {w: i for i, w in enumerate(cols)}
coord = []
for j, idx in enumerate(S):
    for w, c in rows[idx][1]:
        ci = ids[w]
        k = bisect.bisect_left(pp, ci)
        if k < len(pp) and pp[k] == ci:
            coord.append((k, j, int(c)))
mc.log(f'B nnz={len(coord)} s={time.time()-t0:.1f}')
# nonsingularity mod p1
Bn = flint.nmod_mat(5042, 5042, mc.P1)
for i, j, v in coord: Bn[i, j] = v % mc.P1
assert Bn.det() != 0
mc.log(f'B nonsingular mod p1 s={time.time()-t0:.1f}')
B = flint.fmpz_mat(5042, 5042)
for i, j, v in coord: B[i, j] = v
b = flint.fmpz_mat(5042, 2)
for i, c in enumerate(pp):
    w = cols[c]
    gr = Grat.get(w)
    if gr: b[i, 0] = int(gr*Lrat)
    gs = Gsq.get(w)
    if gs: b[i, 1] = int(gs*Lsq)
mc.log(f'fmpz B and RHS built s={time.time()-t0:.1f}; starting exact solve')
t1 = time.time()
Y = B.solve(b)
mc.log(f'solve done s={time.time()-t1:.1f}')
out = {}
for comp, (name, L) in enumerate((('rat', Lrat), ('sq', Lsq))):
    y = [F(int(Y[j, comp].numerator), int(Y[j, comp].denominator))/L for j in range(5042)]
    out[name] = y
    nz = sum(1 for a in y if a)
    h = max((max(abs(a.numerator), a.denominator) for a in y), default=1)
    mc.log(f'y_{name}: nonzero={nz} maxheight_digits={len(str(h))}')
with open(str(mc.HERE/'stage2_y.pkl'), 'wb') as f:
    pickle.dump(out, f)

# EXACT verification on ALL 5184 columns
srows = [dict(rows[j][1]) for j in S]
for name in ('rat', 'sq'):
    acc = {}
    y = out[name]
    for j in range(5042):
        a = y[j]
        if not a: continue
        for w, c in srows[j].items():
            acc[w] = acc.get(w, F(0)) + a*c
    Gx = Grat if name == 'rat' else Gsq
    resid = {}
    for w in set(acc) | set(Gx):
        dd = acc.get(w, F(0)) - Gx.get(w, F(0))
        if dd: resid[w] = dd
    mc.log(f'EXACT verify {name}: accumulated_cols={len(acc)} residual_nonzero={len(resid)} s={time.time()-t0:.1f}')
    assert not resid, (name, next(iter(resid.items())))
mc.log('EXACT VERIFICATION PASSED: G_rat and G_sq are in rowspace_Q(T0), witness rows S (all 5184 columns)')
