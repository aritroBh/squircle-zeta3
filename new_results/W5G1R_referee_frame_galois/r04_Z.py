"""r04: claim 4. (i) Recompute Z independently from W5C2's certificates (exact y), the lift1 rows, and the hand
decomposition of r01 (u.b0 = rho + 4e (u=3), -4e (u=5); u.b1 = rho - 4e (u=3), +4e (u=5)).
(ii) Rebuild every generator G_i named in W5G1's Z_certificate.pkl from provenance (B4 rows from r03, verified against
make_relations; A3h4 frame rows regenerated via frame_ds with own domain checks) and verify sum t_i G_i == Z exactly on
all 5184 weight-4 words, and also on random projections mod two primes.
(iii) compare the 19 frame lifts with those of Theorem W4R's certificate."""
from rcommon import *
import pickle, json, random, collections
sys.set_int_max_str_digits(0)
L = 'r04.log'
W5C1 = T + '/W5C1_certificate_scoping'
W5 = words(5); W4 = words(4)
SEC = {'pp': (1, 1), 'pm': (1, -1), 'mp': (-1, 1), 'mm': (-1, -1)}
def chi(g, s): return (s[0] if g in (1, 3) else 1) * (s[1] if g in (2, 3) else 1)
ECO = {(0, 3): 4, (0, 5): -4, (1, 3): -4, (1, 5): 4}
B = pickle.load(open(HERE + '/r03_bases.pkl', 'rb'))
lift1 = pickle.load(open(W5C1 + '/rels_lift1.pkl', 'rb'))
ex = pickle.load(open(HERE + '/r02_extract.pkl', 'rb'))['exact']
Z = {'rat': {}, 'sq': {}}
vals = []
for tag, fam, idx, lam, yr, ys in ex:
    if fam != 'lift1': continue
    lab, r = lift1[idx]
    _, j, w = lab
    r = {W5[i]: F(c) for i, c in r}
    bw = shv(B[1][j], {tuple(w): 1})
    w0 = next(iter(r)); kap = r[w0] / F(bw[w0])
    assert all(r.get(x, 0) == kap * bw.get(x, 0) for x in set(r) | set(bw))
    for g, u in ((2, 3), (3, 5)):
        for comp, y in (('rat', yr), ('sq', ys)):
            if not y: continue
            a = y / (4 * lam) * chi(g, SEC[tag])
            uw = act(tuple(w), u)
            Z[comp][uw] = Z[comp].get(uw, 0) + a * kap * ECO[(j, u)]
for c in Z: Z[c] = {w: v for w, v in Z[c].items() if v}
log(f'(i) Z recomputed from certificates: support rat {len(Z["rat"])}, sq {len(Z["sq"])} (W5G1 s05: 662 / 666)', L)
# (ii)
ZC = pickle.load(open(T + '/W5G1_frame_galois/Z_certificate.pkl', 'rb'))
sys.path.insert(0, T + '/A3h4_frame_ds'); sys.path.insert(0, T + '/A3c_claude_elementary/Q')
import frame_ds as fd
from cmzv import POINTS, cross
from qsqrt2 import QS2 as Q
from functools import reduce
from itertools import product
def qf(x): return float(x.p) + float(x.q) * 2 ** .5
def regen_frame(n, i):
    data = pickle.load(open(T + f'/A3h4_frame_ds/validated_w{n}.pickle', 'rb'))
    g, fs = data['provenance'][i]
    p0, p1, pinf = g
    ts = {cross(p, p0, p1, pinf): p for p in POINTS}
    bs = {1 / t for t in ts if t is not None and abs(qf(t)) > 1}
    # own domain checks: every factor argument real with 0<|z|<1; every cumulative product and every
    # vertex product of cumulative products is 1/t for a frame point t with |t|>1
    zs = [z for f in fs for k, z in f]
    assert all(0 < abs(qf(z)) < 1 for z in zs)
    cums = []
    for f in fs:
        b = Q(1); s = []
        for k, z in f: b = b * z; s.append(b)
        assert all(x in bs for x in s)
        cums.append(s)
    for vertex in product(*((Q(1),) + tuple(s) for s in cums)):
        if any(x != Q(1) for x in vertex): assert reduce(lambda a, b: a * b, vertex, Q(1)) in bs
    raw = fd.primitive(fd.convert(g, ts, fd.identity(fs)))
    assert raw == data['rows'][i]
    return {tuple(w): F(c) for w, c in raw}, max(abs(qf(z)) for z in zs)
FRC = {}
for comp in ('rat', 'sq'):
    t = ZC[comp]
    acc = {}
    nfr = collections.Counter(); maxz = 0
    for labk, tv in t.items():
        if labk[0] == 'A3h5_S': G = B[4][labk[1]]
        else:
            _, n, i, w = labk
            if (n, i) not in FRC: FRC[(n, i)] = regen_frame(n, i)
            G = shv(FRC[(n, i)][0], {tuple(w): 1}); nfr[(n, i, tuple(w))] += 1; maxz = max(maxz, FRC[(n, i)][1])
        for w, c in G.items(): acc[w] = acc.get(w, 0) + tv * c
    acc = {w: v for w, v in acc.items() if v}
    mism = sum(1 for w in set(acc) | set(Z[comp]) if acc.get(w, 0) != Z[comp].get(w, 0))
    log(f'(ii) {comp}: {len(t)} generators ({sum(nfr.values())} frame lifts, max |z| {maxz:.3f}); EXACT sum t_i G_i == Z on all words: mismatches {mism}', L)
    # random projections mod p
    for p in (2305843009213693951, 1000000007):
        rnd = random.Random(p); pi = {w: rnd.randrange(p) for w in W4}
        mp_ = lambda fr: fr.numerator % p * pow(fr.denominator % p, -1, p) % p
        lhs = sum(mp_(tv) * sum(mp_(F(c)) * pi[w] for w, c in (B[4][k[1]] if k[0] == 'A3h5_S' else shv(FRC[(k[1], k[2])][0], {tuple(k[3]): 1})).items()) for k, tv in t.items()) % p
        rhs = sum(mp_(v) * pi[w] for w, v in Z[comp].items()) % p
        log(f'(ii) {comp}: random projection mod {p}: {"equal" if lhs == rhs else "DIFFERENT"}', L)
    frs = sorted(nfr)
    log(f'(ii) {comp}: frame lifts used {frs}', L)
cert = json.load(open(T + '/A3h5_exact_certificate/D_certificate.json'))
w4r = sorted((int(l[1]), int(l[2]), tuple(l[3])) for l in cert['labels'] if l[0] == 'frame' and int(l[1]) in (2, 3))
log(f'(iii) W4R frame lifts of weight-2/3 rows: {w4r}', L)
log(f'(iii) Z frame lifts subset of W4R frame lifts: {set(frs) <= set(w4r)}', L)
pickle.dump(Z, open(HERE + '/r04_Z.pkl', 'wb'))
