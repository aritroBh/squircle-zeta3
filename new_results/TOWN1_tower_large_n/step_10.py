import os
base = '<project>'
p = os.path.join(base, 'research/08_theory/TOW1_tower_law')
with open(os.path.join(p, 'result.md')) as f:
    print(f.read()[:12000])
