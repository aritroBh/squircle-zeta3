"""s05: split W5C2's unfolded identity  D = sum_s sum_i y_i/(4 lam_i) sum_g chi_s(g) g.r_i  into
  (a) literal rows (g=1) and their complex conjugates (u=7)                     -- justified (B),
  (b) u=3,5 images of std rows: literal rows of the std family                  -- justified,
  (c) u=3,5 images of lift rows  = kappa*(rho sh u.w) + kappa*e sh (y sh u.w):   rho-part = literal lift rows,
      e-part collected into  e sh Z,  Z a weight-4 vector,
  (d) u=3,5 images of frame-w5 rows: literal if present in validated_w5.pkl, else collected into X.
Checks EXACTLY that (a)+(b)+(c)+(d) == D on all 46656 words, both components, and that every claimed
literal row really is in the files.  Writes split.pkl with Z, X, the frame (row,u) list and the literal-term list."""
from common import *
from math import gcd
import json, collections
t0 = time.time()
W5 = words(5); IDS = {w: i for i, w in enumerate(W5)}
used = pickle.load(open(HERE + '/used_rows.pkl', 'rb'))
FAMS = ('std', 'lift1', 'lift2old', 'lift3old', 'lift4old')
R = {f: pickle.load(open(W5C1 + f'/rels_{f}.pkl', 'rb')) for f in FAMS}
VW5 = pickle.load(open(W5C1 + '/validated_w5.pkl', 'rb'))
FR = [{tuple(w): c for w, c in row} for row in VW5['rows']]
B = bases(); DEC = pickle.load(open(HERE + '/lift_decomp.pkl', 'rb'))
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
    return tuple((w, a // g) for w, a in items), F(g, den)     # d = mult * key
# literal index: key -> (fam, idx, mult)
LIT = {}
for f in FAMS:
    for i in range(len(R[f])):
        k, m = normkey(rowdict(f, i)); LIT.setdefault(k, (f, i, m))
for i in range(len(FR)):
    k, m = normkey(FR[i]); LIT.setdefault(k, ('framew5', i, m))
LAB = {f: {R[f][i][0]: i for i in range(len(R[f]))} for f in FAMS if f != 'std'}
log(f'literal index built: {len(LIT)} keys ({time.time()-t0:.0f}s)', 's05.log')
# aggregate coefficients per unique row and g
coef = collections.defaultdict(lambda: {c: [F(0)] * 4 for c in ('rat', 'sq')})
for tag, L in used.items():
    sec = SEC[tag]
    for f, i, lam, yr, ys in L:
        for c, y in (('rat', F(yr)), ('sq', F(ys))):
            if not y: continue
            b = y / (4 * F(lam))
            for gi in range(4): coef[(f, i)][c][gi] += b * chi(gi, sec)
log(f'unique used rows: {len(coef)}', 's05.log')
terms = collections.defaultdict(lambda: {'rat': F(0), 'sq': F(0)})   # (fam, idx, u) -> coefficient; u in {1,7}
Z = {'rat': {}, 'sq': {}}; X = {'rat': {}, 'sq': {}}
frame_pairs = {}    # (idx, u) -> 'literal' / 'X'
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
                f2, i2, m2 = LIT[k]
                s = m / m2
                terms[(f2, i2, 1)]['rat'] += cr * s; terms[(f2, i2, 1)]['sq'] += cs * s
                stats[('literal', f, u)] += 1
                if f == 'framew5': frame_pairs[(i, u)] = 'literal'
                continue
            assert f == 'framew5', ('std image not literal', i, u)
            frame_pairs[(i, u)] = 'X'; stats[('X', f, u)] += 1
            addv(X['rat'], gr, cr); addv(X['sq'], gr, cs); continue
        # lift row: label (fam, j, w)
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
        for wy, cy in y.items():          # e sh (y sh u.w): collect y sh u.w into Z
            addv(Z['rat'], mul({wy: 1}, {gw: 1}), cr * kappa * cy)
            addv(Z['sq'], mul({wy: 1}, {gw: 1}), cs * kappa * cy)
        stats[('lift_decomp', f, u, bool(y))] += 1
log('term stats: ' + str(dict(stats)), 's05.log')
log(f'literal/conjugate terms: {len(terms)}; frame (row,u) pairs: ' + str(collections.Counter(frame_pairs.values())), 's05.log')
# exact reassembly
D = {c: {tuple(int(x) for x in k.split(',')): F(v) for k, v in json.load(open(W5C1 + f'/D_{c}.json')).items()} for c in ('rat', 'sq')}
for c in ('rat', 'sq'):
    acc = {}
    for (f, i, u), cc in terms.items():
        if cc[c]: addv(acc, actv(rowdict(f, i), u), cc[c])
    eZ = mul(E1, Z[c]); addv(acc, eZ, 1); addv(acc, X[c], 1)
    bad = sum(1 for w in set(acc) | set(D[c]) if acc.get(w, F(0)) != D[c].get(w, F(0)))
    log(f'EXACT reassembly {c}: D == literal + conj + e sh Z + X ? mismatches={bad}; |supp Z|={sum(1 for v in Z[c].values() if v)} |supp X|={sum(1 for v in X[c].values() if v)}', 's05.log')
    assert bad == 0
pickle.dump(dict(terms={k: dict(v) for k, v in terms.items()}, Z=Z, X=X, frame_pairs=frame_pairs), open(HERE + '/split.pkl', 'wb'))
log(f'saved split.pkl ({time.time()-t0:.0f}s)', 's05.log')
