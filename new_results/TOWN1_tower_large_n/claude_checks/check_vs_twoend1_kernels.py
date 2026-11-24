# M1 session: compare TOWN1's degree-2 closed forms and degree-3 kernel polynomials with the INDEPENDENT full-reducer
# kernels of TWOEND1 at k = 19, 20, 21 (TWOEND1_weight_poles/red/kern_K.pkl, built by fastred from the full integrand),
# plus WTS's kern_K.pkl at k = 16..18. Exact Fractions.
import pickle, json
from fractions import Fraction as F
R = '<project>/research/08_theory'
law = json.load(open(f'{R}/TOWN1_tower_large_n/results/kernel_laws.json'))['degree3_residue_polynomials']
POLE = {'-a': (F(-1), F(-1)), 'a': (F(1), F(1)), '0': (F(0), F(0)), '-b': (F(1), F(-1)), 'b': (F(-1), F(1))}
def closed2(k):
    # {(key): {pole: (rational, sqrt2)}}
    return {(2, 0): {POLE['0']: (F(1, 4), F(-1, 4)), POLE['a']: (F(1, 4), F(0)), POLE['-a']: (F(2*k - 7, 4), F(0)), POLE['-b']: (F(-(k - 2), 2), F(0))},
            (1, 1): {POLE['0']: (F(2*k - 5, 4), F(1, 4)), POLE['-a']: (F(-(2*k - 5), 2), F(0))}}
def poly3(k):
    out = {}
    for key, kk in [('30', (3, 0)), ('21', (2, 1))]:
        out[kk] = {}
        for pole, d in law[key].items():
            r = sum(F(c)*k**n for n, c in enumerate(d['rat'])); s = sum(F(c)*k**n for n, c in enumerate(d['sqrt2']))
            if r or s: out[kk][POLE[pole]] = (r, s)
    return out
def real_part(P):
    return {p: (c[0], c[1]) for p, c in P['real'].items() if c[0] or c[1]}, all(not any(c) for c in P['real'].values() for _ in [0] if False)
bad = 0; n = 0
srcs = [(k, f'{R}/WTS_weight_k_structure/theory/kern_{k}.pkl') for k in (16, 17, 18)] + [(k, f'{R}/TWOEND1_weight_poles/red/kern_{k}.pkl') for k in (19, 20, 21)]
for k, fn in srcs:
    K = pickle.load(open(fn, 'rb'))
    want = {**closed2(k), **poly3(k)}
    for key, W in want.items():
        P = K[key]
        if isinstance(P, dict) and 'real' in P:
            got = {p: (c[0], c[1]) for p, c in P['real'].items() if c[0] or c[1]}
            imag = any(c[2] or c[3] for c in P['real'].values()) or any(any(e) for e in P['ipos']) or any(any(e) for e in P['ineg'])
        else:
            got = {(r.p0, r.p1): (c.p0, c.p1) for r, c in P[0].items() if c.p0 or c.p1}
            imag = any(c.p2 or c.p3 for c in P[0].values())
        ok = (got == W) and not imag
        n += 1; bad += (not ok)
        print(f'k={k} key={key}: {"MATCH" if ok else "MISMATCH"}' + ('' if ok else f'  got={got} want={W} imag={imag}'))
print(f'{n - bad}/{n} match', 'ALL MATCH' if not bad else 'MISMATCHES')
