import os
tow1 = '<project>/research/08_theory/TOW1_tower_law'
with open(os.path.join(tow1, 'proof_tower.md')) as f:
    content = f.read()
print(len(content))
print(content[:8000])
