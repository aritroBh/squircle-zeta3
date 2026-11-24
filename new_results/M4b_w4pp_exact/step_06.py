import os
ROOT = '<project>'
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
M4A = os.path.join(ROOT,'research/08_theory/M4a_decomposition')
# head of gpl_w4pp.tsv
with open(os.path.join(M4A,'gpl_w4pp.tsv')) as f:
    lines = f.read().splitlines()
print(len(lines))
for l in lines[:15]:
    print(repr(l))
print('...')
import json
print(json.load(open(os.path.join(M4A,'quad_90.json')))['W4pp'])
