import os
ROOT = '<project>'
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
PS = os.path.join(ROOT,'paper/squircle_note/scripts')
src = open(os.path.join(PS,'w4_symbol.py')).read()
print(src)
