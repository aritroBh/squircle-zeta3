import os
base = '<project>/research/08_theory'
tow3 = base + '/TOW3_tower_depth6'
with open(tow3 + '/proof_depth6.md') as f:
    txt = f.read()
print(len(txt))
print(txt[:9000])
