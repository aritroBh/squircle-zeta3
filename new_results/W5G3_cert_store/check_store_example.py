"""Independent spot checker written from FORMAT.md only (no W5G2/W5G3 scripts imported).
Usage: check_store_example.py u3.f4   (reads store/<name>.core.jsonl.xz and store/<name>.lit.jsonl.xz)   -> prints mismatch counts of identity (I5) and (I4)."""
import sys, lzma, json, pickle, itertools
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
R = '<project>/research/08_theory/'
sys.path[:0] = [R + 'A3h4_frame_ds', '<project>/paper/squircle_note/scripts']   # only so pickles of sqrt2 numbers load
W5C1 = R + 'W5C1_certificate_scoping'
def words(n): return [w for w in itertools.product(range(-1, 8), repeat=n) if w[-1] != -1 and w[0] != 0]
W5, W4 = words(5), words(4)
def act(w, u): return tuple(-1 if x == -1 else (u * x) % 8 for x in w)
def q(s): a, b = s.split('/'); return F(int(a), int(b))
sections = {'rows': {}}
for part in ('core', 'lit'):
    with lzma.open(f'store/{sys.argv[1]}.{part}.jsonl.xz', 'rt') as f:
        for line in f:
            o = json.loads(line)
            if o['section'] == 'rows': sections['rows'][o['family']] = o['terms']
            elif o['section'] == 'meta' and 'meta' in sections: pass
            else: sections[o['section']] = o
meta = sections['meta']; i, u = meta['i'], meta['u']
V = pickle.load(open(W5C1 + '/validated_w5.pkl', 'rb'))['rows']
acc = {}
def addto(d, w, c):
    d[w] = d.get(w, 0) + c
# target v
v = {}
for w, c in V[i]: addto(v, act(tuple(w), u), c)
# frame terms
F234 = None
for t in sections['frame']['terms']:
    a = q(t[-1])
    if t[0] == 'framew5':
        for w, c in V[t[1]]: addto(acc, tuple(w), a * c)
    else:
        if F234 is None: F234 = pickle.load(open(W5C1 + '/rels_frame234.pkl', 'rb'))
        for k, c in F234[t[1]][1]: addto(acc, W5[k], a * c)
# literal rows
for fam, terms in sections['rows'].items():
    RL = pickle.load(open(W5C1 + f'/rels_{fam}.pkl', 'rb'))
    for idx, g, s in terms:
        a = q(s)
        for k, c in RL[idx][1]: addto(acc, act(W5[k], g), a * c)
# e sh W
def shuffle_insert(letter, w):
    return [w[:p] + (letter,) + w[p:] for p in range(len(w) + 1)]
for w4, s in sections['W']['terms']:
    a = q(s)
    for x in shuffle_insert(2, tuple(w4)): addto(acc, x, a)
    for x in shuffle_insert(6, tuple(w4)): addto(acc, x, -a)
bad5 = sum(1 for w in set(acc) | set(v) if acc.get(w, 0) != v.get(w, 0))
print(meta['name'], 'I5 mismatches on', len(W5), 'words (support checked on union):', bad5)
# I4: W = sum t_k G_k
Wd = {tuple(w): q(s) for w, s in sections['W']['terms']}
T0 = pickle.load(open(R + 'A3h_weight4_target/w4_systems_T.pkl', 'rb'))['T0']
S = pickle.load(open(R + 'A3h5_exact_certificate/stage2_rows_S.pkl', 'rb'))['S']
def sh(u_, v_):
    if not u_: return {v_: 1}
    if not v_: return {u_: 1}
    out = {}
    for w, a in sh(u_[1:], v_).items(): out[(u_[0],) + w] = out.get((u_[0],) + w, 0) + a
    for w, a in sh(u_, v_[1:]).items(): out[(v_[0],) + w] = out.get((v_[0],) + w, 0) + a
    return out
acc4 = {}
for t in sections['t']['terms']:
    a = q(t[-1])
    if t[0] == 'A3h5_S':
        for w, (qq, c) in T0[S[t[1]]][1]:
            assert qq == 0; addto(acc4, tuple(w), a * F(c))
    else:
        n, ri, wd = t[1], t[2], tuple(t[3])
        row = pickle.load(open(R + f'A3h4_frame_ds/validated_w{n}.pickle', 'rb'))['rows'][ri]
        for w, c in row:
            for x, m in sh(tuple(w), wd).items(): addto(acc4, x, a * c * m)
bad4 = sum(1 for w in set(acc4) | set(Wd) if acc4.get(w, 0) != Wd.get(w, 0))
print(meta['name'], 'I4 mismatches:', bad4)
