import os
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)

with open(os.path.join(base, 'research/08_theory/WTS_weight_k_structure/theory/result.md')) as f:
    wts = f.read()
print(wts[:14000])
