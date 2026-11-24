import os
base = '<project>'
p = os.path.join(base, 'research/08_theory/WT3_weight3_all_k')
with open(os.path.join(p, 'result.md')) as f:
    wt3 = f.read()
print(wt3[:11000])
