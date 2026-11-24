"""t06: every frame relation used by the certificates (t04 frame terms, t05 weight-4 frame lifts) is regenerated from
its stored provenance with A3h4's frame_ds (unchanged), with the domain conditions re-asserted exactly as in
Theorem W4R's verifier (A3h5 verify_D.py): the factor tuple must be produced by candidates(n, bs, arity) of its
frame (which enforces real arguments b = 1/t, |t| > 1, decreasing |b|, every vertex product a frame point), and
primitive(convert(identity(fs))) must equal the stored row.  Also checks u=7 images used are of such rows."""
from g2common import *
import frame_ds as fd
from cmzv import POINTS, cross
LOG = 't06.log'
t0 = time.time()
cs = pickle.load(open(HERE + '/cert_summary.pkl', 'rb'))['summary']
need5 = set(); need234 = set()
f234 = pickle.load(open(W5C1 + '/rels_frame234.pkl', 'rb'))
for rec in cs:
    for (fam, i, u), a in rec['frame_terms']:
        if fam == 'framew5': need5.add(i)
        else:
            lab = f234[i][0]; assert lab[0] == 'frame'; need234.add((lab[1], lab[2]))
            # the lift row is primitive(base sh w): check
            base = pickle.load(open(A3H4 + f'/validated_w{lab[1]}.pickle', 'rb'))['rows'][lab[2]]
            W = words(5); IDS = {w: k for k, w in enumerate(W)}
            row = primitive(mul(dict(base), {tuple(lab[3]): 1}))
            assert {IDS[w]: c for w, c in row} == dict(f234[i][1]), ('frame234 lift mismatch', i)
if os.path.exists(HERE + '/w4cert_summary.pkl'):
    for (n, i) in pickle.load(open(HERE + '/w4cert_summary.pkl', 'rb'))['frame_rows_used']: need234.add((n, i))
log(f'frame-w5 rows used: {sorted(need5)}', LOG)
log(f'A3h4 frame rows (weight, index) used via lifts: {sorted(need234)}', LOG)
def regen(n, g, fs, stored):
    p0, p1, pinf = g
    ts = {cross(p, p0, p1, pinf): p for p in POINTS}
    bs = {1 / t for t in ts if t is not None and abs(t) > 1}
    assert fs in tuple(fd.candidates(n, bs, len(fs))), 'factors not a domain candidate'
    raw = primitive(fd.convert(g, ts, fd.identity(fs)))
    return raw == primitive(dict(stored)), max(abs(float(z.p) + float(z.q) * 2 ** .5) for f in fs for k, z in f)
VW5 = load_framew5(); ok5 = []
for i in sorted(need5):
    g, fs = VW5['provenance'][i]
    ok, mz = regen(5, g, fs, [(tuple(w), c) for w, c in VW5['rows'][i]])
    ok5.append(ok); log(f'  frame-w5 row {i}: regenerated in domain, equals stored: {ok}; max|z| {mz:.4f}', LOG)
ok4 = []
for n, i in sorted(need234):
    data = pickle.load(open(A3H4 + f'/validated_w{n}.pickle', 'rb'))
    g, fs = data['provenance'][i]
    ok, mz = regen(n, g, fs, data['rows'][i])
    ok4.append(ok); log(f'  A3h4 frame row (w{n}, {i}): regenerated in domain, equals stored: {ok}; max|z| {mz:.4f}', LOG)
log(f'ALL frame rows regenerated from provenance within the domain: {all(ok5) and all(ok4)} '
    f'({len(ok5)} weight-5, {len(ok4)} weight 2-4) ({time.time()-t0:.0f}s)', LOG)
assert all(ok5) and all(ok4)
