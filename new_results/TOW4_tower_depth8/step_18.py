import os
base = '<project>/research/08_theory'
with open(base + '/TOW2_tower_all_depths/proof_all_depths.md') as f:
    txt = f.read()
print(len(txt))
# find section 5
idx = txt.find('Section 5')
if idx == -1: idx = txt.find('## 5')
if idx == -1: idx = txt.find('master')
print(txt[max(0,idx-200):idx+6500])
