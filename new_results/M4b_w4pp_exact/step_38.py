import os, subprocess, time
Z3 = 'Zeta' + '(3)'
ROOT = os.path.join('<home>/Downloads/Research', Z3)
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
PY = os.path.join(ROOT, '.venv/bin/python')
src = open('step_03.py').read()   # my pipeline file was written correctly by the earlier (failed) block?
print('len of my step_03.py:', len(src))
print(src[:80])
