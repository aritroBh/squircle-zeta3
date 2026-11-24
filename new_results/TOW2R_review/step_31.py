import os
wts = '<project>/research/08_theory/WTS_weight_k_structure/theory'
with open(os.path.join(wts, 'result.md')) as f:
    content = f.read()
print(content[13300:16200])
