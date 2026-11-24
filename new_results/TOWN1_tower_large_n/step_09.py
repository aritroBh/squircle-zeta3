import os
base = '<project>'
p = os.path.join(base, 'research/08_theory/WT2_weight2_all_k')
for fn in ['buildintegrand.py', 'genintegrand.py']:
    fp = os.path.join(p, fn)
    if os.path.exists(fp):
        with open(fp) as f:
            print("="*30, fn)
            print(f.read())
