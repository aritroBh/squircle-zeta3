import os, subprocess, time
Z3 = 'Zeta' + '(3)'
ROOT = os.path.join('<home>/Downloads/Research', Z3)
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
PY = os.path.join(ROOT, '.venv/bin/python')
step07 = r'''"""step07: write W4pp_certificate.json.
Certificate: W4pp_X = sum_k c_k^X B_k + sum_j z_j^X F_j + sum_i y_i^X T0_{S_i}, X in {rat, sq},
as exact rational word-vector identities (5184 columns), proven by stage1 (mod-p CRT, hold-out)
and stage2 (exact fmpz solve + full Fraction verification).
labels: ['atom', k] (36), ['frame', n, i, w] (28), ['T0', lab] (5040 used rows of S).
coefficients: [rat, sq] as strings, aligned with labels."""
import sys, json, pickle, time
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/M4b_w4pp_exact')
import m4b_common as mc
from fractions import Fraction as F

t0 = time.time()
AT = pickle.load(open(str(mc.HERE/'stage2_atoms.pkl'), 'rb'))['atoms']
d = json.load(open(str(mc.HERE/'stage1_exact_cz.json')))
cr = [F(x) for x in d['c_rat']]; cs = [F(x) for x in d['c_sq']]
zr = [F(x) for x in d['z_rat']]; zs = [F(x) for x in d['z_sq']]
chosen = d['chosen_frames']
frames = mc.load_frame_rows()
Y = pickle.load(open(str(mc.HERE/'stage2_y.pkl'), 'rb')))
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
         "Stage2: exact fmpz solve, full 5184-column Fraction verification (M4b, A3h5 template).",
 'target': 'gpl_w4pp.tsv (M4a_decomposition); word vectors built by step01_build_vectors.py',
 'atoms': [enc for _, enc, _ in AT],
 'atom_names': [n for n, _, _ in AT],
 'labels': labels,
 'coefficients': [car, casq],
}
(mc.HERE/'W4pp_certificate.json').write_text(json.dumps(cert))
from collections import Counter
fam = Counter(l[0] for l in labels)
mc.log(f'W4pp_certificate.json written: labels={len(labels)} families={dict(fam)} '
       f'nonzero_rat={sum(1 for a in car if F(a))} nonzero_sq={sum(1 for a in casq if F(a))} '
       f'bytes={(mc.HERE/"W4pp_certificate.json").stat().st_size} s={time.time()-t0:.1f}')
'''
open('step07_make_certificate.py','w').write(step07)
r = subprocess.run([PY, 'step07_make_certificate.py'], capture_output=True, text=True, cwd=WORK)
print('EXIT', r.returncode)
print(r.stdout)
print(r.stderr[-2000:] if r.returncode else '')
