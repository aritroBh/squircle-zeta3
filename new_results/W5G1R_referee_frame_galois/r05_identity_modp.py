"""r05: claims 1 (usage) and 6. Own code, mod two primes, random projections pi on the 46656 words.
 (a) W5C2's unfolded identity D = sum_entries c_entry sum_g chi(g) g.r, for c = y/(4 lam) and for c = y lam/4.
 (b) X: u=3,5 images of the frame-w5 rows that are not literally frame-w5 rows (exact, own computation).
 (c) D - X - e sh Z  ==  sum of PROVED vectors only:  literal rows (u=1), conjugates (u=7), std images (literal std rows,
     r03), lift images rebuilt as kappa * sum_i t_i (b_i sh u.w) from the proved bases (r03-verified decompositions),
     and literal frame images.  Checked on projections mod p."""
from rcommon import *
import pickle, json, random, collections
from math import gcd, lcm
L = 'r05.log'
W5C1 = T + '/W5C1_certificate_scoping'
sys.path.insert(0, ROOT + '/paper/squircle_note/scripts')
W5 = words(5); IDX = {w: i for i, w in enumerate(W5)}
SEC = {'pp': (1, 1), 'pm': (1, -1), 'mp': (-1, 1), 'mm': (-1, -1)}
UN = (1, 7, 3, 5)
def chi(g, s): return (s[0] if g in (1, 3) else 1) * (s[1] if g in (2, 3) else 1)
P = (2305843009213693951, 1000000007)
FAMS = ('std', 'lift1', 'lift2old', 'lift3old', 'lift4old')
R = {f: pickle.load(open(W5C1 + f'/rels_{f}.pkl', 'rb')) for f in FAMS}
VW5 = pickle.load(open(W5C1 + '/validated_w5.pkl', 'rb'))
FR = [{tuple(w): F(c) for w, c in row} for row in VW5['rows']]
def rowd(f, i): return FR[i] if f == 'framew5' else {W5[k]: F(c) for k, c in R[f][i][1]}
ext = pickle.load(open(HERE + '/r02_extract.pkl', 'rb'))
B = pickle.load(open(HERE + '/r03_bases.pkl', 'rb'))
DEC = pickle.load(open(T + '/W5G1_frame_galois/lift_decomp.pkl', 'rb'))
Z = pickle.load(open(HERE + '/r04_Z.pkl', 'rb'))
D = {c: {tuple(int(x) for x in k.split(',')): F(v) for k, v in json.load(open(W5C1 + f'/D_{c}.json')).items()} for c in ('rat', 'sq')}
def prim(d):
    d = {w: F(c) for w, c in d.items() if c}
    den = lcm(*[c.denominator for c in d.values()])
    it = sorted((w, int(c * den)) for w, c in d.items())
    g = 0
    for _, a in it: g = gcd(g, abs(a))
    if it[0][1] < 0: g = -g
    return tuple((w, a // g) for w, a in it)
FRSET = {prim(r): k for k, r in enumerate(FR)}
# ---- (b) X exact
coefx = collections.defaultdict(lambda: {'rat': F(0), 'sq': F(0)})
for tag, fam, idx, lam, yr, ys in ext['exact']:
    if fam != 'framew5': continue
    for g in (2, 3):
        coefx[(idx, UN[g])]['rat'] += yr / (4 * lam) * chi(g, SEC[tag])
        coefx[(idx, UN[g])]['sq'] += ys / (4 * lam) * chi(g, SEC[tag])
X = {'rat': {}, 'sq': {}}; lit = {}; nonlit = []
for (idx, u), cc in sorted(coefx.items()):
    im = actv(FR[idx], u)
    k = FRSET.get(prim(im))
    if k is not None: lit[(idx, u)] = k; continue
    nonlit.append((idx, u))
    for c in ('rat', 'sq'):
        for w, v in im.items(): X[c][w] = X[c].get(w, 0) + cc[c] * v
for c in X: X[c] = {w: v for w, v in X[c].items() if v}
log(f'(b) frame (row,u) pairs u=3,5: literal {sorted(lit.items())}; non-literal {len(nonlit)} over rows {sorted({i for i, _ in nonlit})}', L)
log(f'(b) X support rat {len(X["rat"])}, sq {len(X["sq"])} (W5G1 s05: 26286 / 26264)', L)
lit7 = sum(1 for i in {i for i, _ in coefx} if prim(actv(FR[i], 7)) in FRSET)
log(f'(b) frame-w5 rows used whose u=7 image is literally a frame-w5 row: {lit7}/{len({i for i, _ in coefx})}', L)
E = {(2,): 1, (6,): -1}
for p in P:
    rnd = random.Random(12345 + p % 1000)
    piv = [rnd.randrange(p) for _ in W5]
    perm = {u: [IDX[act(w, u)] for w in W5] for u in UN}
    def m(fr): fr = F(fr); return fr.numerator % p * pow(fr.denominator % p, -1, p) % p
    def proj(d): return sum(m(c) * piv[IDX[w]] for w, c in d.items()) % p
    rowproj = {}
    def rp(f, i, u):
        k = (f, i, u)
        if k not in rowproj:
            if f == 'framew5': rowproj[k] = proj(actv(FR[i], u))
            else: rowproj[k] = sum(int(c) * piv[perm[u][j]] for j, c in R[f][i][1]) % p
        return rowproj[k]
    for c, ci in (('rat', 4), ('sq', 5)):
        Dp = proj(D[c])
        for conv in ('y/(4lam)', 'y*lam/4'):
            s = 0
            for tag, fam, idx, lam, ymr, yms in ext['mod']:
                y = (ymr if c == 'rat' else yms)[P.index(p)]
                if not y: continue
                lm = m(lam)
                a = y * pow(4 * lm, -1, p) if conv == 'y/(4lam)' else y * lm * pow(4, -1, p)
                for g in range(4): s += a * chi(g, SEC[tag]) * rp(fam, idx, UN[g])
            log(f'(a) p={p} {c}: D.pi == unfolded sum with c={conv}: {Dp == s % p}', L)
        # (c) proved-only reassembly
        s = 0
        for tag, fam, idx, lam, ymr, yms in ext['mod']:
            y = (ymr if c == 'rat' else yms)[P.index(p)]
            if not y: continue
            a = y * pow(4 * m(lam), -1, p)
            for g in range(4):
                u = UN[g]; ag = a * chi(g, SEC[tag])
                if u in (1, 7) or fam == 'std': s += ag * rp(fam, idx, u); continue   # u=1,7 all; std literal (r03)
                if fam == 'framew5':
                    if (idx, u) in lit: s += ag * rp(fam, idx, u)
                    continue                                                       # non-literal -> X
                n = {'lift1': 1, 'lift2old': 2, 'lift3old': 3, 'lift4old': 4}[fam]
                lab, rr = R[fam][idx]; j, w = lab[1], tuple(lab[2])
                r = rowd(fam, idx); bw = shv(B[n][j], {w: 1}); w0 = next(iter(r)); kap = r[w0] / F(bw[w0])
                t, yv = DEC[(n, j, u)]
                uw = act(w, u)
                for i2, ti in t.items():
                    s += ag * m(kap * ti) * proj(shv(B[n][i2], {uw: 1}))
            s %= p
        rhs = (proj(D[c]) - proj(X[c]) - proj(shv(E, Z[c]))) % p
        log(f'(c) p={p} {c}: (D - X - e sh Z).pi == (sum of proved vectors).pi : {rhs == s}', L)
pickle.dump(dict(X=X, lit=lit, nonlit=nonlit), open(HERE + '/r05_X.pkl', 'wb'))
