import os, subprocess
ROOT = '<project>'
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
# find w4_symbol.py and m3_pslq.py
for root, dirs, files in os.walk(ROOT):
    for f in files:
        if f in ('w4_symbol.py','m3_pslq.py'):
            print(os.path.join(root,f))
