import os, json
Z3 = 'Zeta' + '(3)'
ROOT = os.path.join('<home>/Downloads/Research', Z3)
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
d = json.load(open('stage1_exact_cz.json'))
from fractions import Fraction as F
print('selected:', d['selected'])
print('atom names:', len(d['atom_names']))
names = d['atom_names']
cr = [F(x) for x in d['c_rat']]; cs = [F(x) for x in d['c_sq']]
zr = [F(x) for x in d['z_rat']]; zs = [F(x) for x in d['z_sq']]
print('== atom coefficients (rat, sq) ==')
for k,(n,a,b) in enumerate(zip(names, cr, cs)):
    if a or b:
        print(f'{k:3d} {n:28s} {a} , {b}')
print('== frame z (nonzero) ==')
for k,(a,b) in enumerate(zip(zr, zs)):
    if a or b:
        print(k, d['frame_labels'][k], a, b)
