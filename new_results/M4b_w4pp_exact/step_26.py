import os, json
ROOT = '<project>'
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
M4A = os.path.join(ROOT,'research/08_theory/M4a_decomposition')
om = json.load(open(os.path.join(M4A,'omega4.json')))
for k,v in om.items():
    print(k, v)
