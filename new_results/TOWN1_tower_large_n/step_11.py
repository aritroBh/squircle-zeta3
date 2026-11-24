import os
base = '<project>'
p = os.path.join(base, 'research/08_theory/TOW1_tower_law')
with open(os.path.join(p, 'proof_tower.md')) as f:
    proof = f.read()
print(len(proof))
print(proof[:13000])
