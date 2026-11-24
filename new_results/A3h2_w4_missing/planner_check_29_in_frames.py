"""Planner check (A3h2 step 5): how many of the 29 missing weight-4 relations are absorbed by span(T3 + A3h4 frame lifts)?
Uses A3h's mod-p RREF of T3 (sqrt2 -> r) and the frame lifts built as in A3h4_frame_ds/planner_check_membership.py."""
import sys, json, pickle
from pathlib import Path
HERE = Path(__file__).resolve().parent; THEORY = HERE.parent
sys.path.insert(0, str(THEORY / 'A3c_claude_elementary/Q'))
from cmzv import mul
from relations import words, primitive
from fractions import Fraction as F
from flint import nmod_mat

cols = words(4); ids = {w: i for i, w in enumerate(cols)}
J = json.load(open(HERE / 'w4_missing_relations.json'))
fw = [tuple(w) for w in J['freewords']]
rels = J['relations']
print('relations:', len(rels), '; free words:', len(fw))
frame = []
for n in (2, 3, 4):
    for row in pickle.load(open(THEORY / f'A3h4_frame_ds/validated_w{n}.pickle', 'rb'))['rows']:
        for w in words(4 - n):
            lifted = primitive(mul(dict(row), {w: 1}))
            if lifted: frame.append(lifted)


def fm(c, p):
    c = F(c); return c.numerator % p * pow(c.denominator % p, -1, p) % p


for p in (1000000009, 1000000033):
    R = pickle.load(open(THEORY / f'A3h_weight4_target/rref_T3_{p}.pkl', 'rb'))
    pp, free, rowsF = R['pp'], R['free'], R['rowsF']
    fpos = {c: k for k, c in enumerate(free)}; piv = {c: i for i, c in enumerate(pp)}
    def nf(pairs):
        out = [0] * len(free)
        for w, c in pairs:
            a = fm(c, p)
            if not a: continue
            col = ids[tuple(w)]
            if col in fpos: out[fpos[col]] = (out[fpos[col]] + a) % p
            else:
                for f, r in rowsF[piv[col]]: out[fpos[f]] = (out[fpos[f]] - a * r) % p
        return out
    def rank(vs):
        vs = [v for v in vs if any(v)]
        return nmod_mat(len(vs), len(free), [x for v in vs for x in v], p).rank() if vs else 0
    nF = [nf(r) for r in frame]
    nR = [nf([(tuple(int(x) for x in k.split('.')), c) for k, c in r.items()]) for r in rels]
    rF = rank(nF); rR = rank(nR); rFR = rank(nF + nR)
    inside = sum(1 for v in nR if rank(nF + [v]) == rF)
    print(f'p={p}: T3 quotient={len(free)}; rank frames={rF}; rank 29 rels={rR}; rank frames+rels={rFR} '
          f'(=> {rF + rR - rFR}-dim overlap); relations individually inside span(frames): {inside}/29; '
          f'quotient after frames+29 = {len(free) - rFR}')
