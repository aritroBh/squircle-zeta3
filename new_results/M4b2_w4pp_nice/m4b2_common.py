"""Shared loaders for the M4b2 W4'' NICE exact certificate (adapted from A3h5 common.py / planner_probe_z.py).
No results trusted from pickles without independent recomputation; loaders only bring in data."""
import sys, pickle, time, json
sys.dont_write_bytecode = True
from pathlib import Path
HERE = Path('<project>/research/08_theory/M4b2_w4pp_nice')
THEORY = HERE.parent
Z3 = THEORY.parent.parent
for p in (str(Z3/'paper/squircle_note/scripts'), str(THEORY/'A3_W4_proof'),
          str(THEORY/'A3c_claude_elementary/Q'), str(THEORY/'A3h4_frame_ds')):
    if p not in sys.path: sys.path.insert(0, p)
from fractions import Fraction as F
from collections import Counter

P1, P2 = 1000000009, 1000000033
PRIMES = [1000000009, 1000000033, 1000000097, 1000000241, 1000000289, 1000000297, 1000000321, 1000000409]

def log(s):
    s = time.strftime('%Y-%m-%d %H:%M:%S ') + s
    print(s, flush=True)
    with (HERE/'progress.log').open('a') as f: f.write(s + '\n')

def load(path):
    with Path(path).open('rb') as f: return pickle.load(f)

def words4():
    from relations import words
    return words(4)

def load_T0():
    T = load(THEORY/'A3h_weight4_target/w4_systems_T.pkl')['T0']
    rows = []
    for lab, row in T:
        assert all(q == 0 for w, (q, c) in row)
        rows.append((('T0', lab), tuple((w, F(c)) for w, (q, c) in row)))
    return rows

def load_frame_rows():
    """Lifted frame rows in deterministic order n=(2,3,4), label ('frame',n,i,w)."""
    from relations import words, primitive
    from cmzv import mul
    out = []
    for n in (2, 3, 4):
        data = load(THEORY/f'A3h4_frame_ds/validated_w{n}.pickle')
        for i, row in enumerate(data['rows']):
            for w in words(4-n):
                lifted = primitive(mul(dict(row), {w: 1}))
                if lifted:
                    out.append((('frame', n, i, w), tuple((t, F(c)) for t, c in lifted)))
    return out

def fmod(c, p):
    c = F(c)
    return c.numerator % p * pow(c.denominator % p, -1, p) % p

def rref_path(p):
    a = THEORY/f'A3h5_exact_certificate/rref_T0_{p}.pkl'
    b = THEORY/f'A3h_weight4_target/rref_T0_{p}.pkl'
    if a.exists(): return a
    assert b.exists()
    return b

def nf_at(p):
    """Normal form in the T0 quotient at prime p (planner_probe_z style). Returns (nf, m)."""
    R = load(rref_path(p)); assert R['rank'] == 5042 and R['p'] == p
    free, pp, rowsF = R['free'], R['pp'], R['rowsF']
    fpos = {c: k for k, c in enumerate(free)}; piv = {c: i for i, c in enumerate(pp)}
    ids = {w: i for i, w in enumerate(words4())}
    def nf(pairs):
        out = [0] * len(free)
        items = pairs.items() if isinstance(pairs, dict) else pairs
        for w, c in items:
            a = fmod(c, p)
            if not a: continue
            col = ids[w]
            if col in fpos: out[fpos[col]] = (out[fpos[col]] + a) % p
            else:
                for f, r in rowsF[piv[col]]: out[fpos[f]] = (out[fpos[f]] - a * r) % p
        return out
    return nf, len(free)
