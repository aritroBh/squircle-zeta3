import os, json
ROOT = '<project>'
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
M4A = os.path.join(ROOT,'research/08_theory/M4a_decomposition')
print(open(os.path.join(M4A,'gexpr.py')).read())
print('=====CHECK JSON=====')
print(open(os.path.join(M4A,'gpl_w4pp_check.json')).read())
