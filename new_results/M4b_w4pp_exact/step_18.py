import os
ROOT = '<project>'
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
src = open(os.path.join(ROOT,'paper/squircle_note/scripts/qsqrt2.py')).read()
print(src[:8000])
