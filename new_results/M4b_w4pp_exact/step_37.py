import os
Z3 = 'Zeta' + '(3)'
ROOT = os.path.join('<home>/Downloads/Research', Z3)
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
print('workaround ok:', ROOT)
import glob
f = sorted(glob.glob('step_*.py'), key=os.path.getmtime)[-1]
print(repr(open(f).read()[:120]))
