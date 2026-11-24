import os
base = '<project>/research/08_theory'
tow3 = base + '/TOW3_tower_depth6'
with open(tow3 + '/proof_depth6.md') as f:
    txt = f.read()
# print chars 4500..6500 (the truncated section 3)
print(txt[4300:6800])
