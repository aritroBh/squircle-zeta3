"""r02a: regenerate every lift row (lift1, lift2old, lift3old, lift4old) with my own shuffle from the base rows
B1..B4 (data), and every frame234 lift used by W5G2's frame basis from A3h4 validated rows; compare with the
stored rows used by W5G2 (W5C1 rels_*.pkl). Save my rows for r02."""
from rc import *
LOG = 'r02a.log'
t0 = time.time()
B = bases()
mine = {}
for fam, n in (('lift1', 1), ('lift2old', 2), ('lift3old', 3), ('lift4old', 4)):
    stored = load_rows(fam)
    st = {lab: d for lab, d in stored}
    out = []; bad = 0; missing = 0
    for j, b in enumerate(B[n]):
        for w in WORDS(5 - n):
            r = prim(shv(b, {w: 1}))
            if not r: continue
            lab = (fam, j, w)
            out.append((lab, r))
            if lab not in st: missing += 1
            elif st[lab] != r: bad += 1
    log(f'{fam}: regenerated {len(out)} rows (stored {len(stored)}); label missing in stored: {missing}; '
        f'mismatching rows: {bad} ({time.time()-t0:.0f}s)', LOG)
    assert len(out) == len(stored) and bad == 0 and missing == 0
    mine[fam] = out
# frame234 lifts used
cs = pickle.load(open(T + '/W5G2_frame_exact/cert_summary.pkl', 'rb'))['summary']
f234 = pickle.load(open(T + '/W5C1_certificate_scoping/rels_frame234.pkl', 'rb'))
W5 = WORDS(5)
used = sorted({lab[1] for rec in cs for lab, a in rec['frame_terms'] if lab[0] == 'frame234'})
A3 = {n: pickle.load(open(T + f'/A3h4_frame_ds/validated_w{n}.pickle', 'rb'))['rows'] for n in (2, 3, 4)}
fr = {}; bad = 0
for i in used:
    lab, row = f234[i]
    _, n, k, w = lab
    r = prim(shv({tuple(x): c for x, c in A3[n][k]}, {tuple(w): 1}))
    if r != {W5[a]: c for a, c in row}: bad += 1
    fr[i] = (lab, r)
log(f'frame234 lifts in the frame bases: {len(used)} distinct; regenerated with own shuffle; mismatches {bad}; '
    f'base rows {sorted({(l[1], l[2]) for l, _ in fr.values()})}', LOG)
assert bad == 0
mine['frame234'] = fr
pickle.dump(mine, open(HERE + '/my_lifts.pkl', 'wb'))
log(f'saved my_lifts.pkl ({time.time()-t0:.0f}s)', LOG)
