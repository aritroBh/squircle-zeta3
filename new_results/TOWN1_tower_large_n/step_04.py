import os
base = '<project>'
with open(os.path.join(base, 'research/08_theory/WTS_weight_k_structure/theory/result.md')) as f:
    wts = f.read()
print(wts[14000:26000])
