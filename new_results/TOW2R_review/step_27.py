import os
base = '<project>/research/08_theory'
wts = os.path.join(base, 'WTS_weight_k_structure', 'theory')
files = sorted(os.listdir(wts))
print([f for f in files if f.startswith('tower_')][:30])
print([f for f in files if f.startswith('kern_')])
# WTS result.md
with open(os.path.join(wts, 'result.md')) as f:
    content = f.read()
print(len(content))
