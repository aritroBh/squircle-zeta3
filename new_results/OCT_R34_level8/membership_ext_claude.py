"""R3/R4 DIFF vectors (from verify_oct_r34_claude.py) vs base + Kummer (544) and + R1..R5 (549), mod two primes.
Also: rank of base + Kummer + DIFF_Lam + DIFF_swap, and which R_i combinations the DIFFs produce.
Run from research/08_theory after verify_oct_r34_claude.py."""
import sys, pickle
from fractions import Fraction as Fr
import flint
from sympy.ntheory import sqrt_mod
sys.argv = ['x']
exec(open('OCT_R34_level8/verify_oct_r34_claude.py').read().split("Wrows = [accum(r)")[0])   # atoms, targets, loaders
Wrows = [accum(r) for r in W['rows']]
DIFFS = {}
for key, (forms, _) in targets.items():
    G = C.forms_path(forms, P[Q(0)], P[B]); R = expr_vec(red[key])
    DIFFS[key] = {w: c for w, c in C.add((G, 1), (R, -1)).items() if c != 0}
F5 = pickle.load(open('A3f4_glanois_basis/five_relations.pickle', 'rb'))
free, rels = F5['freewords'], F5['relations']
for p in (1000000033, 1000000241):
    s = int(sqrt_mod(2, p))
    base = []
    for r in Wrows:
        v = [0] * n
        for w, a in r.items(): v[idx[w]] = (v[idx[w]] + modp(a, p)) % p
        base.append(v)
    for r in K:
        v = [0] * n
        for w, (x, y) in r.items(): v[idx[w]] = (v[idx[w]] + modp(x, p) + s * modp(y, p)) % p
        base.append(v)
    Rv = []
    for rel in rels:
        v = [0] * n
        for w, a in zip(free, rel): v[idx[tuple(w)]] = (v[idx[tuple(w)]] + modp(a, p)) % p
        Rv.append(v)
    dv = {k: [modp(D.get(c, 0), p) for c in cols] for k, D in DIFFS.items()}
    rk = lambda rows: flint.nmod_mat(len(rows), n, [x for r in rows for x in r], p).rank()
    b = rk(base); bR = rk(base + Rv)
    out = [f"p={p}: base+K {b}, +R1..R5 {bR}"]
    for k, d in dv.items():
        out.append(f"{k}: +DIFF {rk(base + [d])}, +R+DIFF {rk(base + Rv + [d])}")
    out.append(f"both DIFFs: {rk(base + list(dv.values()))}")
    for i in range(5):
        out.append(f"base+both DIFFs+R{i+1}: {rk(base + list(dv.values()) + [Rv[i]])}")
    print(" | ".join(out), flush=True)
