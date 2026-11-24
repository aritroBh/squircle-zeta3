"""Planner check of A3h4's key claim, using A3h's (Kimi's) mod-p RREF of T0, not Codex's elimination code.

Claim: with Codex's numerically validated frame double-shuffle rows (weights 2, 3, 4, lifted to weight 4 by shuffling
with every word of the complementary weight), D_rat and D_sq lie in span(T0 + frame lifts) + Q*P4 mod p.

nf(v)[f] = v[f] - sum_i v[pp_i] * R_i[f] (mod p), where row i of the RREF is e_{pp_i} + sum_f R_i[f] e_f.
Run from this folder:  PY -B planner_check_membership.py
"""
import sys, pickle
from pathlib import Path
HERE = Path(__file__).resolve().parent
THEORY = HERE.parent
sys.path.insert(0, str(THEORY / 'A3c_claude_elementary/Q'))
from cmzv import mul
from relations import words, primitive
from fractions import Fraction as F
from flint import nmod_mat

A3H = THEORY / 'A3h_weight4_target'
cols = words(4); ids = {w: i for i, w in enumerate(cols)}
T = pickle.load(open(A3H / 'w4_systems_T.pkl', 'rb'))['T0']
V = pickle.load(open(A3H / 'w4_systems.pkl', 'rb'))

frame = []
for n in (2, 3, 4):
    rows = pickle.load(open(HERE / f'validated_w{n}.pickle', 'rb'))['rows']
    for row in rows:
        for w in words(4 - n):
            lifted = primitive(mul(dict(row), {w: 1}))
            if lifted: frame.append(lifted)
print('frame lifted rows:', len(frame))


def fm(c, p):
    c = F(c); return c.numerator % p * pow(c.denominator % p, -1, p) % p


for p in (1000000009, 1000000033):
    R = pickle.load(open(A3H / f'rref_T0_{p}.pkl', 'rb'))
    pp, free, rowsF = R['pp'], R['free'], R['rowsF']
    fpos = {c: k for k, c in enumerate(free)}
    piv = {c: i for i, c in enumerate(pp)}

    def nf(vec):  # vec: dict col -> int mod p
        out = [0] * len(free)
        for c, a in vec.items():
            if c in fpos: out[fpos[c]] = (out[fpos[c]] + a) % p
            else:
                for f, r in rowsF[piv[c]]: out[fpos[f]] = (out[fpos[f]] - a * r) % p
        return out

    def vec(pairs): return {ids[w]: fm(c, p) for w, c in pairs if fm(c, p)}
    # sanity: T0 rows reduce to 0
    bad = sum(1 for lab, row in T[:400] if any(nf(vec([(w, cc) for w, (q, cc) in row]))))
    nfF = [nf(vec(r)) for r in frame]
    nP = nf(vec([(w, c) for w, (q, c) in V['P4'].items()]))
    nA = nf(vec([(w, c) for w, (q, c) in V['Drat'].items()]))
    nB = nf(vec([(w, c) for w, (q, c) in V['Dsq'].items()]))

    def rank(vs):
        vs = [v for v in vs if any(v)]
        return nmod_mat(len(vs), len(free), [x for v in vs for x in v], p).rank() if vs else 0
    rF = rank(nfF); rFP = rank(nfF + [nP])
    print(f'p={p}: T0 sanity bad={bad}/400; rank nf(frames)={rF}; +P4={rFP}; +P4+Drat={rank(nfF+[nP,nA])}; '
          f'+P4+Dsq={rank(nfF+[nP,nB])}; +P4+Drat+Dsq={rank(nfF+[nP,nA,nB])}; T0 quotient={len(free)} -> {len(free)-rF}')


# ---- coefficient of P4 mod p: reduce nf(D) and nf(P4) modulo span(nf(frames)); then residual(D) = c * residual(P4)
def reduce_mod_span(basis_rows, v, p):
    v = v[:]
    for piv_c, row in basis_rows:
        a = v[piv_c]
        if a:
            for j, x in enumerate(row):
                if x: v[j] = (v[j] - a * x) % p
    return v


for p in (1000000009, 1000000033):
    R = pickle.load(open(A3H / f'rref_T0_{p}.pkl', 'rb'))
    pp, free, rowsF = R['pp'], R['free'], R['rowsF']
    fpos = {c: k for k, c in enumerate(free)}; piv = {c: i for i, c in enumerate(pp)}
    def nf(vec):
        out = [0] * len(free)
        for c, a in vec.items():
            if c in fpos: out[fpos[c]] = (out[fpos[c]] + a) % p
            else:
                for f, r in rowsF[piv[c]]: out[fpos[f]] = (out[fpos[f]] - a * r) % p
        return out
    def vec(pairs): return {ids[w]: fm(c, p) for w, c in pairs if fm(c, p)}
    vs = [nf(vec(r)) for r in frame]
    M = nmod_mat(len(vs), len(free), [x for v in vs for x in v], p)
    E, rk = M.rref()
    basis = []
    for i in range(rk):
        row = [int(E[i, j]) for j in range(len(free))]
        pc = next(j for j, x in enumerate(row) if x)
        basis.append((pc, row))
    rP = reduce_mod_span(basis, nf(vec([(w, c) for w, (q, c) in V['P4'].items()])), p)
    out = []
    for name in ('Drat', 'Dsq'):
        rD = reduce_mod_span(basis, nf(vec([(w, c) for w, (q, c) in V[name].items()])), p)
        j = next(k for k, x in enumerate(rP) if x)
        c = rD[j] * pow(rP[j], -1, p) % p
        resid = sum(1 for x, y in zip(rD, rP) if (x - c * y) % p)
        out.append((name, c, resid))
    print(f'p={p}: P4 coefficient mod p (component, c, leftover nonzeros):', out)
