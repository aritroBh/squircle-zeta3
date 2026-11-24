"""WT4 cert step07 (usage: wt4_cert_step07.py <D1|D2>); adapted from M4b2 step07 (NICE): write W4pp_nice_certificate.json.
Certificate: W4pp_X = sum_k c_k^X B_k + sum_j z_j^X F_j + sum_i y_i^X T0_{S_i}, X in {rat, sq},
as exact rational word-vector identities (5184 columns), proven by stage1 (mod-p CRT, hold-out)
and stage2 (exact fmpz solve + full Fraction verification).
labels: ['atom', k] (36), ['frame', n, i, w] (28), ['T0', lab] (5040 used rows of S).
coefficients: [rat, sq] as strings, aligned with labels."""
import sys, json, pickle, time
sys.dont_write_bytecode = True
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))
import wt4_cert_common as mc
from fractions import Fraction as F

t0 = time.time()
TGT = sys.argv[1]
assert TGT in ('D1','D2')
AT = pickle.load(open(os.path.join(HERE, 'stage2_atoms_nice.pkl'), 'rb'))['atoms']
d = json.load(open(os.path.join(HERE, f'stage1_exact_cz_{TGT}.json')))
cr = [F(x) for x in d['c_rat']]; cs = [F(x) for x in d['c_sq']]
zr = [F(x) for x in d['z_rat']]; zs = [F(x) for x in d['z_sq']]
chosen = d['chosen_frames']
frames = mc.load_frame_rows()
Y = pickle.load(open(os.path.join(HERE, f'stage2_y_{TGT}.pkl'), 'rb'))
rows = mc.load_T0()
S = pickle.load(open(str(mc.THEORY/'A3h5_exact_certificate/stage2_rows_S.pkl'), 'rb'))['S']

labels = []; car = []; casq = []
for k, (name, enc, v) in enumerate(AT):
    if not (cr[k] or cs[k]): continue
    labels.append(['atom', k]); car.append(str(cr[k])); casq.append(str(cs[k]))
for j, idx in enumerate(chosen):
    if not (zr[j] or zs[j]): continue
    lab = frames[idx][0]
    labels.append(['frame', lab[1], lab[2], list(lab[3])]); car.append(str(zr[j])); casq.append(str(zs[j]))
used = [j for j in range(len(S)) if Y['rat'][j] or Y['sq'][j]]
for j in used:
    lab = rows[S[j]][0]
    labels.append(['T0', [lab[1][0], lab[1][1]] + ([list(lab[1][2])] if len(lab[1]) > 2 else [])])
    car.append(str(Y['rat'][j])); casq.append(str(Y['sq'][j]))

cert = {
 'note': "W4''_X = sum coeffs[X] . rows, X in {rat=0, sq=1}; W4'' = R'' + sqrt2-part; "
         "rows: atom (classical weight-4 constants via cmzv.atom/mul), frame (function-level double shuffle lifts), "
         "T0 (standard weight-4 level-8 relations). Stage1: mod-p CRT at 8 primes with hold-out. "
         "Stage2: exact fmpz solve, full 5184-column Fraction verification (M4b2 (NICE atom set), M4b/A3h5 template).",
 'target': f'{TGT.lower()}.tsv (WT4_weight4_all_k); word vectors built by wt4_cert_step01.py',
 'atoms': [enc for _, enc, _ in AT],
 'atom_names': [n for n, _, _ in AT],
 'labels': labels,
 'coefficients': [car, casq],
}
open(os.path.join(HERE, f'{TGT}_certificate.json'), 'w').write(json.dumps(cert))
from collections import Counter
fam = Counter(l[0] for l in labels)
mc.log(f'[{TGT}] {TGT}_certificate.json written: labels={len(labels)} families={dict(fam)} '
       f'nonzero_rat={sum(1 for a in car if F(a))} nonzero_sq={sum(1 for a in casq if F(a))} '
       f'bytes={os.path.getsize(os.path.join(HERE, f'{TGT}_certificate.json'))} s={time.time()-t0:.1f}')
