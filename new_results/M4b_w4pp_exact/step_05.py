import os
ROOT = '<project>'
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
print(open(os.path.join(ROOT,'research/08_theory/M4a_decomposition/result.md')).read())
