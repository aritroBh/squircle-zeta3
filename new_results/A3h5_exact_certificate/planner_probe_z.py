"""Planner probe: Stage-1 coefficients z (D_X = sum z_j F_j mod T0, for the 28 chosen frame lifts) at every prime
with a saved T0 RREF; then CRT and rational reconstruction, with a hold-out prime as the check. Read-only on data."""
import json
from math import isqrt
from fractions import Fraction as F
from flint import nmod_mat
from common import *

st = json.load(open(HERE / 'stage1_mod.json'))
chosen = st['chosen']
frames = load_frame_rows(); T = load_targets()
Fsel = [frames[i][1] for i in chosen]
primes = [1000000009, 1000000033, 1000000097, 1000000241, 1000000289, 1000000297, 1000000321, 1000000409]


def nf_at(p):
    path = HERE / f'rref_T0_{p}.pkl'
    if not path.exists(): path = THEORY / f'A3h_weight4_target/rref_T0_{p}.pkl'
    R = load(path); assert R['rank'] == 5042
    free, pp, rowsF = R['free'], R['pp'], R['rowsF']
    fpos = {c: k for k, c in enumerate(free)}; piv = {c: i for i, c in enumerate(pp)}
    ids = {w: i for i, w in enumerate(words4())}
    def nf(pairs):
        out = [0] * len(free)
        for w, c in pairs:
            a = fmod(c, p)
            if not a: continue
            col = ids[w]
            if col in fpos: out[fpos[col]] = (out[fpos[col]] + a) % p
            else:
                for f, r in rowsF[piv[col]]: out[fpos[f]] = (out[fpos[f]] - a * r) % p
        return out
    return nf, len(free)


def solve_z(p):
    nf, m = nf_at(p)
    cols = [nf(Fr) for Fr in Fsel]
    res = []
    for name in ('Drat', 'Dsq'):
        b = nf(list(T[name].items()))
        # augmented system M z = b, M is m x 28; solve via rref of [M | b]
        A = nmod_mat(m, 29, [x for i in range(m) for x in [cols[j][i] for j in range(28)] + [b[i]]], p)
        E, rk = A.rref()
        assert rk == 28, (p, name, rk)  # consistent and full column rank
        z = [0] * 28
        for i in range(28):
            row = [int(E[i, j]) for j in range(29)]
            pc = next(j for j in range(28) if row[j])
            assert pc == i
            z[i] = row[28]
        res.append(z)
    return res


def recon(a, m):
    bound = isqrt(m // 2); r0, r1 = m, a; t0, t1 = 0, 1
    while abs(r1) > bound:
        q = r0 // r1; r0, r1 = r1, r0 - q * r1; t0, t1 = t1, t0 - q * t1
    if t1 == 0 or abs(t1) > bound: return None
    return F(r1, t1) if t1 > 0 else F(-r1, -t1)


Z = {p: solve_z(p) for p in primes}
print('z mod p computed at', len(primes), 'primes')
for k, name in enumerate(('Drat', 'Dsq')):
    for use in (len(primes) - 1,):
        ps = primes[:use]; hold = primes[use]
        M = 1; res = [0] * 28
        for p in ps:
            zp = Z[p][k]
            inv = pow(M % p, -1, p) if M > 1 else 1
            res = [r + M * ((zp[i] - r) * inv % p) for i, r in enumerate(res)] if M > 1 else zp[:]
            M *= p
        rec = [recon(r % M, M) for r in res]
        ok = all(c is not None and fmod(c, hold) == Z[hold][k][i] for i, c in enumerate(rec))
        hts = [max(abs(c.numerator), c.denominator) for c in rec if c is not None]
        print(f'{name}: CRT over {use} primes ({M.bit_length()} bits): reconstructed {sum(c is not None for c in rec)}/28, '
              f'hold-out prime agrees: {ok}; max height digits: {max(len(str(h)) for h in hts) if hts else None}')
        if ok: print('   z =', [str(c) for c in rec])

# save the reconstructed z (both components) for the resumed run
out = {'chosen': chosen, 'labels': [repr(frames[i][0]) for i in chosen]}
for k, name in enumerate(('Drat', 'Dsq')):
    M = 1; res = None
    for p in primes:
        zp = Z[p][k]
        if res is None: res = zp[:]; M = p; continue
        inv = pow(M % p, -1, p); res = [r + M * ((zp[i] - r) * inv % p) for i, r in enumerate(res)]; M *= p
    out['z_' + name] = [str(recon(r % M, M)) for r in res]
json.dump(out, open(HERE / 'stage1_exact_z.json', 'w'), indent=1)
print('wrote stage1_exact_z.json')
