import os
k1f = '<project>/research/08_theory/K1f_hexagon_order5'
with open(os.path.join(k1f, 'step_24.py')) as f:
    content = f.read()
# print first 120 lines
lines = content.split('\n')
for i, ln in enumerate(lines[:130]):
    print(f"{i:3d}: {ln}")
