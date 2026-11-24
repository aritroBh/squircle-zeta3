import os
base = '<project>/research/08_theory'
with open(base + '/TOW2_tower_all_depths/result.md') as f:
    txt = f.read()
print(len(txt))
# find master equation
idx = txt.find('master')
print(txt[max(0,idx-3000):idx+4000])
