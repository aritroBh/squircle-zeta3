"""step04: Stage 1 mod p. In the T0 quotient (142 dims), at each of the 8 saved primes:
  solve nf(W4pp_X) = sum c_k nf(B_k) + sum z_j nf(F_j),  X in {rat, sq},
with an independent column set among [frames(28 chosen by A3h5), atoms(36)], frames first,
the SAME columns at all primes (pivot columns at p1, verified full-rank and consistent at all 8).
CRT over 7 primes + rational reconstruction, 8th prime hold-out. Saves stage1_exact_cz.json."""
import sys, json, pickle, time
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/M4b_w4pp_exact')
import m4b_common as mc
from fractions import Fraction as F
from flint import nmod_mat
from math import isqrt

t0 = time.time()
W = pickle.load(open(str(mc.HERE/'stage1_w4pp_vectors.pkl'), 'rb'))
AT = pickle.load(open(str(mc.HERE/'stage2_atoms.pkl'), 'rb'))['atoms']
z0 = json.load(open(str(mc.THEORY/'A3h5_exact_certificate/stage1_exact_z.json')))
chosen = z0['chosen']
frames = mc.load_frame_rows()
for k, idx in enumerate(chosen):
    assert str(frames[idx][0]) == z0['labels'][k], (k, frames[idx][0], z0['labels'][k])
mc.log(f'{len(chosen)} A3h5-chosen frame lifts verified against stage1 labels')
Fsel = [frames[i][1] for i in chosen]
Bvecs = [v for _, _, v in AT]
cols_all = Fsel + Bvecs
nfcols = len(cols_all)
targets = {'rat': W['rat'], 'sq': W['sq']}
primes = mc.PRIMES

# nf of all columns and targets at all primes
NFS = {}
for p in primes:
    nf, m = mc.nf_at(p)
    M = [nf(list(c.items()) if isinstance(c, dict) else list(c)) for c in cols_all]
    T = {x: nf(list(t.items())) for x, t in targets.items()}
    NFS[p] = (M, T, m)
mc.log(f'nf computed at 8 primes, dim={m}, s={time.time()-t0:.1f}')

# independent column set at p1, frames first
p1 = primes[0]
M1 = NFS[p1][0]; m = NFS[p1][2]
A = nmod_mat(m, nfcols, [M1[j][i] for i in range(m) for j in range(nfcols)], p1)
E, rk = A.rref()
sel = []
start = 0
for i in range(rk):
    while not E[i, start]: start += 1
    sel.append(start); start += 1
mc.log(f'column rank at p1: {rk} of {nfcols} ([28 frames | 36 atoms]); selected={sel}')
nsel = len(sel)
assert all(s < 28 or s >= 28 for s in sel)

def solve_at(p):
    Mp, Tp, m = NFS[p]
    nrhs = 2
    A = nmod_mat(m, nsel + nrhs,
                 [x for i in range(m)
                  for x in [Mp[sel[j]][i] for j in range(nsel)] + [Tp['rat'][i], Tp['sq'][i]]], p)
    E, rk = A.rref()
    # pivots must be exactly the first nsel columns
    piv = []; start = 0
    for i in range(rk):
        while not E[i, start]: start += 1
        piv.append(start); start += 1
    assert rk == nsel and piv == list(range(nsel)), (p, rk, piv[:40])
    zr = [int(E[i, nsel]) for i in range(nsel)]
    zs = [int(E[i, nsel+1]) for i in range(nsel)]
    return zr, zs

Z = {}
for p in primes:
    zr, zs = solve_at(p)
    Z[p] = (zr, zs)
mc.log(f'solved at all 8 primes: selected columns independent and systems consistent everywhere')

def recon(a, m):
    bound = isqrt(m // 2); r0, r1 = m, a; t0, t1 = 0, 1
    while abs(r1) > bound:
        q = r0 // r1; r0, r1 = r1, r0 - q*r1; t0, t1 = t1, t0 - q*t1
    if t1 == 0 or abs(t1) > bound: return None
    return F(r1, t1) if t1 > 0 else F(-r1, -t1)

def crt(ps, vals):
    M = 1; res = None
    for p, v in zip(ps, vals):
        if res is None: res = list(v); M = p; continue
        inv = pow(M % p, -1, p)
        res = [r + M*((vi - r)*inv % p) for r, vi in zip(res, v)]
        M *= p
    return res, M

out = {'chosen_frames': chosen, 'frame_labels': [z0['labels'][k] for k in range(len(chosen))],
       'atom_names': [name for name, _, _ in AT], 'selected': sel,
       'primes': primes, 'note': 'W4pp_X = sum_k c_k[X] B_k + sum_j z_j[X] F_j mod T0; selected indexes [28 frames|36 atoms]'}
for comp, k in (('rat', 0), ('sq', 1)):
    ps = primes[:7]; hold = primes[7]
    res, M = crt(ps, [Z[p][k] for p in ps])
    rec = [recon(r % M, M) for r in res]
    ok = all(c is not None and mc.fmod(c, hold) == Z[hold][k][i] for i, c in enumerate(rec))
    hts = [max(abs(c.numerator), c.denominator) for c in rec if c is not None]
    mc.log(f'{comp}: CRT over 7 primes ({M.bit_length()} bits): reconstructed {sum(c is not None for c in rec)}/{nsel}, '
           f'hold-out agrees: {ok}; max height digits={max((len(str(h)) for h in hts), default=0)}')
    assert ok
    # full vector over all 64 columns (0 on unselected)
    full = [F(0)]*nfcols
    for j, s in enumerate(sel): full[s] = rec[j]
    out['c_rat' if comp == 'rat' else 'c_sq'] = [str(c) for c in full[28:]]
    out['z_rat' if comp == 'rat' else 'z_sq'] = [str(c) for c in full[:28]]
json.dump(out, open(str(mc.HERE/'stage1_exact_cz.json'), 'w'), indent=1)
mc.log(f'wrote stage1_exact_cz.json s={time.time()-t0:.1f}')
