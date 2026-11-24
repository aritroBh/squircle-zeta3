"""M4b2 step04: Stage 1 (mod 8 primes, CRT, hold-out) for the NICE atom set.
Columns = [28 A3h5 frame lifts | 29 atoms (22 products + 7 Li4 at the Theorem-W4 arguments)].
Solve W4pp_X = sum c_k^X B_k + sum z_j^X F_j in the T0 quotient at each prime; CRT x7 + hold-out.
Then compare c_k EXACTLY against the M4a coefficients (w4pp_closed_form.json), term by term.
If the atom columns are dependent modulo (frames + T0) the c_k are not unique; then instead
FIX c_k := M4a coefficients and solve only for z_j (28 unknowns), checking consistency at every
prime (rank 28, no pivot into the RHS), CRT x7 + hold-out. Writes stage1_exact_cz_nice.json."""
import sys, json, pickle, time
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/M4b2_w4pp_nice')
import m4b2_common as mc
from fractions import Fraction as F
from math import isqrt
from flint import nmod_mat

t0 = time.time()
W = pickle.load(open(str(mc.HERE/'stage1_w4pp_vectors.pkl'), 'rb'))
AT = pickle.load(open(str(mc.HERE/'stage2_atoms_nice.pkl'), 'rb'))['atoms']
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
NF = 28  # number of frame columns

NFS = {}
for p in primes:
    nf, m = mc.nf_at(p)
    M = [nf(list(c.items()) if isinstance(c, dict) else list(c)) for c in cols_all]
    T = {x: nf(list(t.items())) for x, t in targets.items()}
    NFS[p] = (M, T, m)
mc.log(f'nf computed at 8 primes, dim={m}, s={time.time()-t0:.1f}')

p1 = primes[0]
M1 = NFS[p1][0]; m = NFS[p1][2]
A = nmod_mat(m, nfcols, [M1[j][i] for i in range(m) for j in range(nfcols)], p1)
E, rk = A.rref()
sel = []; start = 0
for i in range(rk):
    while not E[i, start]: start += 1
    sel.append(start); start += 1
mc.log(f'column rank at p1: {rk} of {nfcols} ([28 frames | 29 atoms]); selected={sel}')
unselected_atoms = [s - NF for s in range(NF, nfcols) if s not in sel]
mc.log(f'unselected (dependent) atom indexes: {unselected_atoms} '
       f'({[AT[k][0] for k in unselected_atoms]})')

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

def solve_selected_at(p, sel):
    Mp, Tp, m = NFS[p]
    nsel = len(sel)
    A = nmod_mat(m, nsel + 2,
                 [x for i in range(m)
                  for x in [Mp[sel[j]][i] for j in range(nsel)] + [Tp['rat'][i], Tp['sq'][i]]], p)
    E, rk = A.rref()
    piv = []; start = 0
    for i in range(rk):
        while not E[i, start]: start += 1
        piv.append(start); start += 1
    assert rk == nsel and piv == list(range(nsel)), (p, rk, piv[:40])
    return [int(E[i, nsel]) for i in range(nsel)], [int(E[i, nsel+1]) for i in range(nsel)]

def crt_holdout(sol, nsel):
    out = {}
    for comp, k in (('rat', 0), ('sq', 1)):
        ps = primes[:7]; hold = primes[7]
        res, M = crt(ps, [sol[p][k] for p in ps])
        rec = [recon(r % M, M) for r in res]
        ok = all(c is not None and mc.fmod(c, hold) == sol[hold][k][i] for i, c in enumerate(rec))
        hts = [max(abs(c.numerator), c.denominator) for c in rec if c is not None]
        mc.log(f'  {comp}: CRT over 7 primes ({M.bit_length()} bits): reconstructed '
               f'{sum(c is not None for c in rec)}/{nsel}, hold-out agrees: {ok}; '
               f'max height digits={max((len(str(h)) for h in hts), default=0)}')
        assert ok
        out[comp] = rec
    return out

# ---- M4a reference coefficients, aligned to our atom order ----
m4a = json.load(open(str(mc.THEORY/'M4a_decomposition/w4pp_closed_form.json')))['closed_form_terms']
name2idx = {name: k for k, (name, _, _) in enumerate(AT)}
cm4a = {'rat': [F(0)]*len(AT), 'sq': [F(0)]*len(AT)}
for coeff, expr in m4a:
    if expr.startswith('sqrt2*'):
        comp, name = 'sq', expr[len('sqrt2*'):]
    else:
        comp, name = 'rat', expr
    assert name in name2idx, ('M4a term not an atom', name)
    assert cm4a[comp][name2idx[name]] == 0, ('duplicate', name)
    cm4a[comp][name2idx[name]] = F(coeff)
assert all(any(cm4a[x][k] for x in ('rat', 'sq')) for k in range(len(AT))), 'an atom got no M4a coefficient'
mc.log(f'M4a closed form parsed: {len(m4a)} terms over all {len(AT)} atoms (49 terms expected: {len(m4a)==49})')

out = {'chosen_frames': chosen, 'frame_labels': [z0['labels'][k] for k in range(len(chosen))],
       'atom_names': [name for name, _, _ in AT], 'primes': primes,
       'note': 'W4pp_X = sum_k c_k[X] B_k + sum_j z_j[X] F_j mod T0; NICE atom set (M4b2)'}

comparison = []
def report_comparison(cr, cs, tag):
    nbad = 0
    for k, (name, _, _) in enumerate(AT):
        ok = (cr[k] == cm4a['rat'][k] and cs[k] == cm4a['sq'][k])
        nbad += (not ok)
        comparison.append((name, str(cr[k]), str(cm4a['rat'][k]), str(cs[k]), str(cm4a['sq'][k]), ok))
    mc.log(f'[{tag}] term-by-term comparison vs M4a: {len(AT)*2 - nbad}/{len(AT)*2} components agree; '
           f'mismatches: {nbad}')
    return nbad == 0

method = None
if rk == nfcols:
    mc.log('all 57 columns independent: c_k unique; solving for c_k and z_j (M4b step04 style)')
    sol = {p: solve_selected_at(p, sel) for p in primes}
    rec = crt_holdout(sol, nfcols)
    cr = [rec['rat'][NF + k] for k in range(len(AT))]   # atom part of the unique solution
    cs = [rec['sq'][NF + k] for k in range(len(AT))]
    if report_comparison(cr, cs, 'free 57-column solve'):
        method = 'free_solve_all_independent'
        out['c_rat'] = [str(c) for c in cr]; out['c_sq'] = [str(c) for c in cs]
        # z are the first 28 entries of the full solution (frames are columns 0..27)
        out['z_rat'] = [str(rec['rat'][j]) for j in range(NF)]
        out['z_sq'] = [str(rec['sq'][j]) for j in range(NF)]
    else:
        mc.log('UNEXPECTED: unique solve disagrees with M4a; falling back to fixed-c solve')
else:
    mc.log(f'atom columns DEPENDENT modulo (frames+T0): rank {rk} < {nfcols}; c_k not unique. '
           f'Fixing c_k := M4a coefficients and solving only for z_j.')

if method is None:
    # ---- fixed-c solve: columns = 28 frames; RHS = T - sum_k c_k B_k ----
    cr, cs = cm4a['rat'], cm4a['sq']
    def solve_z_at(p):
        Mp, Tp, m = NFS[p]
        rhs = {}
        for comp in ('rat', 'sq'):
            cc = cr if comp == 'rat' else cs
            v = list(Tp[comp])
            for k in range(len(AT)):
                a = mc.fmod(cc[k], p)
                if not a: continue
                col = Mp[NF + k]
                v = [(vi - a*ci) % p for vi, ci in zip(v, col)]
            rhs[comp] = v
        A = nmod_mat(m, NF + 2,
                     [x for i in range(m)
                      for x in [Mp[j][i] for j in range(NF)] + [rhs['rat'][i], rhs['sq'][i]]], p)
        E, rk2 = A.rref()
        piv = []; start = 0
        for i in range(rk2):
            while not E[i, start]: start += 1
            piv.append(start); start += 1
        # frames independent and system consistent: pivots exactly 0..27, none in RHS columns
        assert rk2 == NF and piv == list(range(NF)), (p, rk2, piv[:32])
        return [int(E[i, NF]) for i in range(NF)], [int(E[i, NF+1]) for i in range(NF)]
    sol = {p: solve_z_at(p) for p in primes}
    mc.log('fixed-c solve: frame columns independent and system CONSISTENT at all 8 primes')
    rec = crt_holdout(sol, NF)
    method = 'fixed_c_solve_z'
    out['c_rat'] = [str(c) for c in cr]; out['c_sq'] = [str(c) for c in cs]
    out['z_rat'] = [str(c) for c in rec['rat']]; out['z_sq'] = [str(c) for c in rec['sq']]

out['method'] = method
out['comparison'] = [{'atom': n, 'c_rat_solved': a, 'c_rat_m4a': b, 'c_sq_solved': c, 'c_sq_m4a': d, 'agree': ok}
                     for (n, a, b, c, d, ok) in comparison]
json.dump(out, open(str(mc.HERE/'stage1_exact_cz_nice.json'), 'w'), indent=1)
mc.log(f'wrote stage1_exact_cz_nice.json (method={method}) s={time.time()-t0:.1f}')
