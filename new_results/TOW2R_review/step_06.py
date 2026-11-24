import os
tow2 = '<project>/research/08_theory/TOW2_tower_all_depths'
with open(os.path.join(tow2, 'tow2lib.py')) as f:
    src = f.read()
# print the truncated part
lines = src.split('\n')
print('\n'.join(lines[100:160]))
