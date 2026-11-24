import os, glob
Z3 = 'Zeta' + '(3)'
ROOT = os.path.join('<home>/Downloads/Research', Z3)
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
for f in ['step_01.py','step_02.py']:
    c = open(f).read()
    print(f, len(c), '->', c.splitlines()[0][:70])
# They look intact (mine). Rename to collision-proof names and verify the pickles they made.
import pickle
W = pickle.load(open('stage1_w4pp_vectors.pkl','rb'))
print('W4pp vectors:', len(W['rat']), len(W['sq']))
AT = pickle.load(open('stage2_atoms.pkl','rb'))
print('atoms:', len(AT['atoms']), 'converts:', len(AT['converts']), 'fails:', AT['fails'])
os.rename('step_01.py','step01_build_vectors.py')
os.rename('step_02.py','step02_build_atoms.py')
print('renamed')
