import os
ROOT = '<project>'
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
CC = os.path.join(ROOT,'research/08_theory/A3f4_glanois_basis/Q/check_conversion.py')
src = open(CC).read()
print(src[:10000])
