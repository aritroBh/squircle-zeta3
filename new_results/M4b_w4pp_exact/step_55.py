import os
Z3 = 'Zeta' + '(3)'
ROOT = os.path.join('<home>/Downloads/Research', Z3)
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
src = open(os.path.join(ROOT,'research/08_theory/A3h4_frame_ds/frame_ds.py')).read()
print(src[:6000])
