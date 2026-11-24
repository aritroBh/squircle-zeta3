"""(1) Exact certificate Z' = sum_i t_i G_i (W5G1 s10's generators: A3h5's 5042 proved weight-4 rows T0[S] and
weight-4 lifts of the A3h4 frame rows of weights 2-4).  Then e sh Z' = sum t_i (e sh G_i) is a relation.
(2) Every frame-w5 (row,u) pair sent to X' in s09 must be one of the 74 pairs certified exactly in W5G2
(cert_summary.pkl), so X' = sum c_(i,u) u.f_i is a combination of proved relations."""
import sys
sys.dont_write_bytecode=True
sys.path.insert(0,'<project>/research/08_theory/W5G1_frame_galois')
from common import *
import common as G1
sys.set_int_max_str_digits(0)
ME='<project>/research/08_theory/W5C14_fourteen_form'
sp = pickle.load(open(ME + '/split_delta.pkl', 'rb'))
# (2) frame pairs
cs = pickle.load(open(T + '/W5G2_frame_exact/cert_summary.pkl', 'rb'))['summary']
certified = {(e['i'], e['u']) for e in cs}
xp = {k for k, v in sp['frame_pairs'].items() if v == 'X'}
lit = {k for k, v in sp['frame_pairs'].items() if v == 'literal'}
log(f'W5G2 certified pairs: {len(certified)}; pairs in X\': {len(xp)}; literal pairs: {sorted(lit)}')
log(f'X\' pairs not certified by W5G2: {sorted(xp - certified)}')
assert xp <= certified
# (1) Z'
B = bases(); W4 = words(4)
gens = list(B[4]); labels = [('A3h5_S', i) for i in range(len(B[4]))]
for n in (2, 3, 4):
    for i, row in enumerate(pickle.load(open(T + f'/A3h4_frame_ds/validated_w{n}.pickle', 'rb'))['rows']):
        for w in words(4 - n):
            gens.append({k: F(v) for k, v in mul(dict(row), {w: 1}).items()}); labels.append(('A3h4_frame', n, i, w))
Zr = {w: c for w, c in sp['Z']['rat'].items() if c}; Zs = {w: c for w, c in sp['Z']['sq'].items() if c}
log(f'|supp Z\'| rat/sq = {len(Zr)}/{len(Zs)}')
xs, rk = exact_combination(gens, [Zr, Zs], W4)
out = {}
for c, x in zip(('rat', 'sq'), xs):
    assert x is not None, c
    t = {labels[i]: v for i, v in enumerate(x) if v}
    out[c] = t
    used_frames = sorted({(k[1], k[2]) for k in t if k[0] == 'A3h4_frame'})
    log(f"Z'_{c} = exact combination of {len(t)} proved weight-4 rows (frame lifts from A3h4 rows {used_frames}); rank of generators {rk}")
pickle.dump(out, open(ME + '/Zprime_certificate.pkl', 'wb'))
log('saved Zprime_certificate.pkl')
