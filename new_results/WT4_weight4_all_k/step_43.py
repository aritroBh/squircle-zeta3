import os
work='<project>/research/08_theory/WT4_weight4_all_k'
step04 = r'''"""WT4 cert step04: Stage 1 (mod 8 primes, CRT, hold-out) for target D1 or D2.
Columns = [28 A3h5 frame lifts | 29 atoms (22 products + 7 Li4 at the Theorem-W4 arguments)].
Solve TGT_X = sum c_k^X B_k + sum z_j^X F_j in the T0 quotient at each prime; CRT x7 + hold-out.
Then compare c_k EXACTLY against the WT4 LLL closed forms (lll_results.json), term by term.
Usage: wt4_cert_step04.py <D1|D2>"""
import sys, json, pickle, time, os
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wt4_cert_common as mc
from fractions import Fraction as F
from math import isqrt
from flint import nmod_mat

t0 = time.time()
TGT = sys.argv[1]
assert TGT in ('D1', 'D2')
W = pickle.load(open(os.path.join(HERE, f'stage1_{TGT}_vectors.pkl'), 'rb'))
AT = pickle.load(open(os.path.join(HERE, 'stage2_atoms_nice.pkl'), 'rb'))['atoms']
z0 = json.load(open(str(mc.THEORY/'A3h5_exact_certificate/stage1_exact_z.json')))
chosen = z0['chosen']
frames = mc.load_frame_rows()
for k, idx in enumerate(chosen):
    assert str(frames[idx][0]) == z0['labels'][k], (k, frames[idx][0], z0['labels'][k])
mc.log(f'[{TGT}] {len(chosen)} A3h5-chosen frame lifts verified against stage1 labels')
Fsel = [frames[i][1] for i in chosen]
Bvecs = [v for _, _, v in AT]
cols_all = Fsel + Bvecs
nfcols = len(cols_all)
targets = {'rat': W['rat'], 'sq': W['sq']}
primes = mc.PRIMES
NF = 28

NFS = {}
for p in primes:
    nf, m = mc.nf_at(p)
    M = [nf(list(c.items()) if isinstance(c, dict) else list(c)) for c in cols_all]
    T = {x: nf(list(t.items())) for x, t in targets.items()}
    NFS[p] = (M, T, m)
mc.log(f'[{TGT}] nf computed at 8 primes, dim={m}, s={time.time()-t0:.1f}')

p1 = primes[0]
M1 = NFS[p1][0]; m = NFS[p1][2]
Amat = nmod_mat(m, nfcols, [M1[j][i] for i in range(m) for j in range(nfcols)], p1)
E, rk = Amat.rref()
sel = []; start = 0
for i in range(rk):
    while not E[i, start]: start += 1
    sel.append(start); start += 1
mc.log(f'[{TGT}] column rank at p1: {rk} of {nfcols} ([28 frames | 29 atoms])')
unselected_atoms = [s - NF for s in range(NF, nfcols) if s not in sel]
mc.log(f'[{TGT}] unselected (dependent) atom indexes: {unselected_atoms} '
       f'({[AT[k][0] for k in unselected_atoms]})')
assert rk == nfcols, 'columns not independent; would need fix-c/solve-z fallback'

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
    Am = nmod_mat(m, nsel + 2,
                 [x for i in range(m)
                  for x in [Mp[sel[j]][i] for j in range(nsel)] + [Tp['rat'][i], Tp['sq'][i]]], p)
    E, rk2 = Am.rref()
    piv = []; start = 0
    for i in range(rk2):
        while not E[i, start]: start += 1
        piv.append(start); start += 1
    assert rk2 == nsel and piv == list(range(nsel)), (p, rk2, piv[:40])
    return [int(E[i, nsel]) for i in range(nsel)], [int(E[i, nsel+1]) for i in range(nsel)]

def crt_holdout(sol, nsel):
    out = {}
    for comp, k in (('rat', 0), ('sq', 1)):
        ps = primes[:7]; hold = primes[7]
        res, M = crt(ps, [sol[p][k] for p in ps])
        rec = [recon(r % M, M) for r in res]
        ok = all(c is not None and mc.fmod(c, hold) == sol[hold][k][i] for i, c in enumerate(rec))
        hts = [max(abs(c.numerator), c.denominator) for c in rec if c is not None]
        mc.log(f'  [{TGT}] {comp}: CRT over 7 primes ({M.bit_length()} bits): reconstructed '
               f'{sum(c is not None for c in rec)}/{nsel}, hold-out agrees: {ok}; '
               f'max height digits={max((len(str(h)) for h in hts), default=0)}')
        assert ok
        out[comp] = rec
    return out

# ---- WT4 LLL reference coefficients, aligned to our atom order ----
lll = json.load(open(os.path.join(HERE, 'lll_results.json')))
cf = lll[f'{TGT}_closed_form']
alias = {"pi^4": "pi^2*pi^2", "pi^2*log2^2": "pi^2*log2*log2", "log2^4": "log2*log2*log2*log2"}
cf = {alias.get(nm, nm): F(v) for nm, v in cf.items()}
name2idx = {name: k for k, (name, _, _) in enumerate(AT)}
cref = {'rat': [F(0)]*len(AT), 'sq': [F(0)]*len(AT)}
for nm, v in cf.items():
    if nm.startswith('sqrt2*'):
        comp, name = 'sq', nm[len('sqrt2*'):]
    else:
        comp, name = 'rat', nm
    assert name in name2idx, ('closed-form term not an atom', name)
    cref[comp][name2idx[name]] = v
mc.log(f'[{TGT}] WT4 LLL closed form parsed: {sum(1 for x in cf.values() if x)} terms over {len(AT)} atoms')

out = {'chosen_frames': chosen, 'frame_labels': [z0['labels'][k] for k in range(len(chosen))],
       'atom_names': [name for name, _, _ in AT], 'primes': primes,
       'note': f'{TGT}_X = sum_k c_k[X] B_k + sum_j z_j[X] F_j mod T0; NICE atom set (WT4, M4b2 template)'}

comparison = []
def report_comparison(cr, cs, tag):
    nbad = 0
    for k, (name, _, _) in enumerate(AT):
        ok = (cr[k] == cref['rat'][k] and cs[k] == cref['sq'][k])
        nbad += (not ok)
        comparison.append((name, str(cr[k]), str(cref['rat'][k]), str(cs[k]), str(cref['sq'][k]), ok))
    mc.log(f'[{TGT}] [{tag}] term-by-term vs WT4 LLL closed form: {len(AT)*2 - nbad}/{len(AT)*2} agree; '
           f'mismatches: {nbad}')
    return nbad == 0

mc.log(f'[{TGT}] all 57 columns independent: c_k unique; solving for c_k and z_j')
sol = {p: solve_selected_at(p, sel) for p in primes}
rec = crt_holdout(sol, nfcols)
cr = [rec['rat'][NF + k] for k in range(len(AT))]
cs = [rec['sq'][NF + k] for k in range(len(AT))]
assert report_comparison(cr, cs, 'free 57-column solve'), f'{TGT}: solved c_k != LLL closed form'
out['c_rat'] = [str(c) for c in cr]; out['c_sq'] = [str(c) for c in cs]
out['z_rat'] = [str(rec['rat'][j]) for j in range(NF)]
out['z_sq'] = [str(rec['sq'][j]) for j in range(NF)]
out['method'] = 'free_solve_all_independent'
out['comparison'] = comparison
json.dump(out, open(os.path.join(HERE, f'stage1_exact_cz_{TGT}.json'), 'w'), indent=1)
mc.log(f'[{TGT}] saved stage1_exact_cz_{TGT}.json s={time.time()-t0:.1f}')
'''
p=os.path.join(work,'wt4_cert_step04.py')
open(p,'w').write(step04)
print('step04 written')
