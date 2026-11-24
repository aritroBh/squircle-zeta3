"""W5G1 s05 applied to Delta: split the unfolded identity
   Delta = sum_s sum_i y_i/(4 lam_i) sum_g chi_s(g) g.r_i
into (a) literal rows and u=7 conjugates, (b) u=3,5 images of std rows (literal std rows),
(c) u=3,5 images of lift rows = literal lift rows + e sh Z' (W5G1's exact per-base-row decompositions lift_decomp.pkl),
(d) u=3,5 images of frame-w5 rows: literal if in validated_w5.pkl, else collected into X'.
Checks EXACTLY that (a)+(b)+(c)+(d) == Delta on all words, both components. Code = W5G1 s05 with inputs swapped."""
import sys
sys.dont_write_bytecode=True
sys.path.insert(0,'<project>/research/08_theory/W5G1_frame_galois')
from common import *          # W5G1 helpers (read-only): words, mul, act, actv, chi, SEC, E1, bases, FAM_N, UNITS
import common as G1
from math import gcd
import json, collections
ME='<project>/research/08_theory/W5C14_fourteen_form'
def lg(s): log(s)            # print only (no file arg -> nothing written to W5G1)
t0 = time.time()
W5 = words(5); IDS = {w: i for i, w in enumerate(W5)}
used = {t: pickle.load(open(ME + f'/used_rows_delta_{t}.pkl', 'rb')) for t in ('pp','pm','mp','mm')}
FAMS = ('std', 'lift1', 'lift2old', 'lift3old', 'lift4old')
R = {f: pickle.load(open(W5C1 + f'/rels_{f}.pkl', 'rb')) for f in FAMS}
VW5 = pickle.load(open(W5C1 + '/validated_w5.pkl', 'rb'))
FR = [{tuple(w): c for w, c in row} for row in VW5['rows']]
B = bases(); DEC = pickle.load(open(G1.HERE + '/lift_decomp.pkl', 'rb'))
def rowdict(f, i):
    if f == 'framew5': return dict(FR[i])
    return {W5[w]: c for w, c in R[f][i][1]}
def normkey(d):
    d = {w: F(c) for w, c in d.items() if c}
    from math import lcm
    den = lcm(*[c.denominator for c in d.values()])
    items = sorted((w, int(c * den)) for w, c in d.items())
    g = 0
    for _, a in items: g = gcd(g, abs(a))
    if items[0][1] < 0: g = -g
    return tuple((w, a // g) for w, a in items), F(g, den)
LIT = {}
for f in FAMS:
    for i in range(len(R[f])):
        k, m = normkey(rowdict(f, i)); LIT.setdefault(k, (f, i, m))
for i in range(len(FR)):
    k, m = normkey(FR[i]); LIT.setdefault(k, ('framew5', i, m))
LAB = {f: {R[f][i][0]: i for i in range(len(R[f]))} for f in FAMS if f != 'std'}
lg(f'literal index built: {len(LIT)} keys ({time.time()-t0:.0f}s)')
coef = collections.defaultdict(lambda: {c: [F(0)] * 4 for c in ('rat', 'sq')})
for tag, L in used.items():
    sec = SEC[tag]
    for f, i, lam, yr, ys in L:
        for c, y in (('rat', F(yr)), ('sq', F(ys))):
            if not y: continue
            b = y / (4 * F(lam))
            for gi in range(4): coef[(f, i)][c][gi] += b * chi(gi, sec)
lg(f'unique used rows: {len(coef)}  families: {collections.Counter(f for f,i in coef)}')
terms = collections.defaultdict(lambda: {'rat': F(0), 'sq': F(0)})
Z = {'rat': {}, 'sq': {}}; X = {'rat': {}, 'sq': {}}
frame_pairs = {}; frame_coef = collections.defaultdict(lambda: {'rat': F(0), 'sq': F(0)})
def addv(acc, d, k):
    for w, c in d.items(): acc[w] = acc.get(w, F(0)) + k * c
stats = collections.Counter()
for (f, i), cf in coef.items():
    r = rowdict(f, i)
    for gi, u in enumerate(UNITS):
        cr, cs = cf['rat'][gi], cf['sq'][gi]
        if not cr and not cs: continue
        if u in (1, 7):
            terms[(f, i, u)]['rat'] += cr; terms[(f, i, u)]['sq'] += cs; stats[('u17', f)] += 1; continue
        gr = actv(r, u)
        if f == 'std' or f == 'framew5':
            k, m = normkey(gr)
            if k in LIT:
                f2, i2, m2 = LIT[k]; s = m / m2
                terms[(f2, i2, 1)]['rat'] += cr * s; terms[(f2, i2, 1)]['sq'] += cs * s
                stats[('literal', f, u)] += 1
                if f == 'framew5': frame_pairs[(i, u)] = 'literal'
                continue
            assert f == 'framew5', ('std image not literal', i, u)
            frame_pairs[(i, u)] = 'X'; stats[('X', f, u)] += 1
            frame_coef[(i, u)]['rat'] += cr; frame_coef[(i, u)]['sq'] += cs
            addv(X['rat'], gr, cr); addv(X['sq'], gr, cs); continue
        n = FAM_N[f]; lab = R[f][i][0]; j, w = lab[1], lab[2]
        bw = mul(B[n][j], {w: 1})
        w0 = next(iter(r)); kappa = F(r[w0]) / F(bw[w0])
        assert all(F(r.get(x, 0)) == kappa * F(bw.get(x, 0)) for x in set(r) | set(bw)), (f, i)
        t, y = DEC[(n, j, u)]
        gw = act(w, u)
        for i2, ti in t.items():
            lab2 = (f, i2, gw)
            if lab2 not in LAB[f]:
                assert not mul(B[n][i2], {gw: 1}), ('missing lift', lab2); continue
            idx2 = LAB[f][lab2]; r2 = rowdict(f, idx2)
            bw2 = mul(B[n][i2], {gw: 1}); x0 = next(iter(r2)); k2 = F(r2[x0]) / F(bw2[x0])
            s = kappa * ti / k2
            terms[(f, idx2, 1)]['rat'] += cr * s; terms[(f, idx2, 1)]['sq'] += cs * s
        for wy, cy in y.items():
            addv(Z['rat'], mul({wy: 1}, {gw: 1}), cr * kappa * cy)
            addv(Z['sq'], mul({wy: 1}, {gw: 1}), cs * kappa * cy)
        stats[('lift_decomp', f, u, bool(y))] += 1
lg('term stats: ' + str(dict(stats)))
lg(f'literal/conjugate terms: {len(terms)}; frame (row,u) pairs: ' + str(collections.Counter(frame_pairs.values())))
lg('frame rows sent to X: ' + str(sorted({i for (i,u),v in frame_pairs.items() if v=='X'})))
dd = pickle.load(open(ME + '/delta.pkl', 'rb'))
D = {c: {w: a for w, a in dd[c].items()} for c in ('rat', 'sq')}
for c in ('rat', 'sq'):
    acc = {}
    for (f, i, u), cc in terms.items():
        if cc[c]: addv(acc, actv(rowdict(f, i), u), cc[c])
    eZ = mul(E1, Z[c]); addv(acc, eZ, 1); addv(acc, X[c], 1)
    bad = sum(1 for w in set(acc) | set(D[c]) if acc.get(w, F(0)) != D[c].get(w, F(0)))
    lg(f'EXACT reassembly {c}: Delta == literal + conj + e sh Z + X ? mismatches={bad}; |supp Z|={sum(1 for v in Z[c].values() if v)} |supp X|={sum(1 for v in X[c].values() if v)}')
    assert bad == 0
pickle.dump(dict(terms={k: dict(v) for k, v in terms.items()}, Z=Z, frame_pairs=frame_pairs,
                 frame_coef={k: dict(v) for k, v in frame_coef.items()}), open(ME + '/split_delta.pkl', 'wb'))
lg(f'saved split_delta.pkl ({time.time()-t0:.0f}s)')
