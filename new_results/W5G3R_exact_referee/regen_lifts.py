"""Regenerate (own shuffle, common.py) every lift row of the families the certificates use and compare exactly:
 - weight-5 lift1, lift2old, lift4old (W5C1 rels_*.pkl): primitive(base row shuffle word);
 - the 32520 frame234 rows: primitive(A3h4 validated row shuffle word);
 - the weight-4 generators A3h5_S = T0[S] (5042 rows): std rows from W4RX's generator (regen_std.pkl),
   lifts with own shuffle.
Bases: B1, B2, B3 = A3c weight{1,2,3}.pickle selections, checked to be among W4RX-regenerated rows
(cumulative, as W4RX check_b); B4 = T0[S] as regenerated here."""
from common import *
sys.path.insert(0, T + '/W4RX_independent_check')
import w4rx_rel as RX
LOG = 'regen_lifts.log'
t0 = time.time()
RS = pickle.load(open(HERE + '/regen_std.pkl', 'rb'))
# lower bases
B = {}; lower = {}
for n in (1, 2, 3):
    rows, fam = RX.standard_rows(n, dict(lower))
    sel = load_data(T + f'/A3c_claude_elementary/Q/weight{n}.pickle')[2]
    sel = [tuple((tuple(w), int(c)) for w, c in r) for r in sel]
    mine = set(rows)
    log(f'B{n}: {len(sel)} selected base rows; all among {len(rows)} W4RX-regenerated weight-{n} rows: {all(r in mine for r in sel)}', LOG)
    assert all(r in mine for r in sel)
    lower[n] = sel; B[n] = [dict(r) for r in sel]
# T0[S]
std4 = RS['w1to4'][4]
TT = load_data(T + '/A3h_weight4_target/w4_systems_T.pkl')['T0']
S = load_data(T + '/A3h5_exact_certificate/stage2_rows_S.pkl')['S']
cnt = {}; bad = []; notprim = 0; B4 = []
for j, i in enumerate(S):
    lab, row = TT[i]
    st = {tuple(w): Fraction(c) for w, (q, c) in row}
    assert all(q == 0 for w, (q, c) in row)
    if lab[0] == 'std': mine = std4[lab[1]]
    else:
        n = {'lift1': 1, 'lift2old': 2, 'lift3old': 3}[lab[0]]
        w = tuple(lab[2]); assert len(w) == 4 - n and w in INDEX(4 - n)
        mine = primitive(shuffle_vec(B[n][lab[1]], {w: 1}))
    cnt[lab[0]] = cnt.get(lab[0], 0) + 1
    if primitive(st) != mine: bad.append(j)
    if tuple(sorted(st.items())) != tuple((w, Fraction(c)) for w, c in mine): notprim += 1
    B4.append(dict(mine))
log(f'T0[S]: {len(S)} rows {cnt}; regenerated (std: W4RX generator; lifts: own shuffle) and equal after primitive normalisation: '
    f'{len(S)-len(bad)}/{len(S)}; stored row literally equal to the primitive regenerated row: {len(S)-notprim}/{len(S)}', LOG)
B[4] = B4
# weight-5 lifts
W5 = WORDS(5)
for fam, n in (('lift1', 1), ('lift2old', 2), ('lift4old', 4)):
    stored = load_data(T + f'/W5C1_certificate_scoping/rels_{fam}.pkl')
    st = {}
    for lab, row in stored:
        assert lab not in st; st[lab] = tuple((W5[k], int(c)) for k, c in row)
    mine = {}
    for j, b in enumerate(B[n]):
        for w in WORDS(5 - n):
            r = primitive(shuffle_vec(b, {w: 1}))
            if r: mine[(fam, j, w)] = r
    same = sum(1 for lab in st if mine.get(lab) == st[lab])
    log(f'{fam}: stored {len(st)} rows, regenerated {len(mine)}; stored rows equal to regenerated row with the same label: {same}/{len(st)}; '
        f'labels only in regenerated: {len(set(mine)-set(st))} ({time.time()-t0:.0f}s)', LOG)
# frame234
A3 = {n: load_data(T + f'/A3h4_frame_ds/validated_w{n}.pickle')['rows'] for n in (2, 3, 4)}
F234 = load_data(T + '/W5C1_certificate_scoping/rels_frame234.pkl')
same = 0
for lab, row in F234:
    _, n, k, w = lab
    assert len(w) == 5 - n and tuple(w) in INDEX(5 - n)
    r = primitive(shuffle_vec({tuple(x): int(c) for x, c in A3[n][k]}, {tuple(w): 1}))
    same += r == tuple((W5[a], int(c)) for a, c in row)
log(f'frame234: {same}/{len(F234)} stored rows equal primitive(A3h4 row shuffle word) with own shuffle ({time.time()-t0:.0f}s)', LOG)
