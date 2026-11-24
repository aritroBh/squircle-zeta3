src = open('gold5_lib.py').read()
old = """def ft_expand_sym_nf(S, dim):
    \"\"\"S: {exp-vector-word -> NF}. Expand letters to basis indices; NF coeffs multiply Fractions.\"\"\"
    out = {}
    for w, cf in S.items():
        terms = [((), Fr(1))]
        for l in w:
            terms = [(p + (j,), c * l[j]) for p, c in terms for j in range(dim) if l[j]]
        for p, c in terms:
            nc = tuple(c * x for x in cf)"""
new = """def ft_expand_sym_nf(S, dim):
    \"\"\"S: {exp-vector-word -> NF or Fr}. Expand letters to basis indices; NF coeffs multiply Fractions.\"\"\"
    out = {}
    for w, cf in S.items():
        if not isinstance(cf, tuple):
            cf = (cf,)
        terms = [((), Fr(1))]
        for l in w:
            terms = [(p + (j,), c * l[j]) for p, c in terms for j in range(dim) if l[j]]
        for p, c in terms:
            nc = tuple(c * x for x in cf)"""
assert old in src
open('gold5_lib.py','w').write(src.replace(old, new))
print("patched")
import sys, itertools, time, json
sys.path.insert(0, '<project>/research/08_theory/GOLD5_golden_corner_w5')
from fractions import Fraction as Fr
import gold5_lib as G
import mpmath as mp
mp.mp.dps = 40
r5 = mp.sqrt(5)
P = G.golden_point('gol1')
F = P.F
ONE, ZERO = F.one(), F.zero()
def fval(x):
    return mp.mpf(x[0].numerator)/x[0].denominator + mp.mpf(x[1].numerator)/x[1].denominator*r5
def vec_or_none(el):
    try: return P.vec(el)
    except ValueError: return None
sols = [ONE, F.neg(ONE)]
for tup in itertools.product(range(-9, 10), repeat=2):
    if tup == (0,0): continue
    v = ONE
    for e_, g in zip(tup, G.BASIS):
        if e_: v = F.mul(v, F.pow(g, e_))
    for sgn in (1, -1):
        x = v if sgn == 1 else F.neg(v)
        if F.eq(x, ONE): continue
        om = F.sub(ONE, x)
        if all(t == 0 for t in om): continue
        if vec_or_none(om) is not None and not any(F.eq(x, y) for y in sols):
            sols.append(x)
atoms = sorted({(x, y) for x in sols for y in sols
                if not F.eq(x, ONE) and any(F.eq(F.mul(x, y), s) for s in sols)
                and abs(fval(x)) <= 1 and abs(fval(F.mul(x, y))) < 1},
               key=lambda t: (float(abs(fval(t[0]))), float(abs(fval(t[1])))))
print(f"{len(atoms)} convergent atoms")
t0 = time.time()
adata = {}
for x, y in atoms:
    dw = G.delta_word_sym(P, G.li32_word(P, x, y), ZERO, ONE)
    dE = G.ft_expand_sym_nf(dw, 2)
    qcw, _ = G.quotient_coords_nf(dE)
    S_at = G.ft_expand_sym_nf(G.li32_sym_nf(P, x, y), 2)
    obs = G.obs_of_SE(S_at)
    adata[(str(x), str(y))] = {'cvec': [[str(c) for c in q] for q in qcw],
                               'obs': [[str(c) for c in o] for o in obs]}
print(f"atom cvecs/obs done ({time.time()-t0:.0f}s)")
json.dump(adata, open('atom_cobrackets_golden.json', 'w'), indent=1)
nz = [(k, v) for k, v in adata.items() if any(Fr(c) for q in v['cvec'] for c in q)]
print(f"{len(nz)} atoms with nonzero cobracket; sample:")
for k, v in nz[:8]:
    print(f"  Li32{k}: cvec={v['cvec']}, obs={v['obs']}")
