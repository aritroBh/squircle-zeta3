import os
base = '<project>'
p = os.path.join(base, 'research/08_theory/TOW1_tower_law')
with open(os.path.join(p, 'towlib.py')) as f:
    towlib = f.read()
# Print the part around fold_feed / Lambda_ij that was truncated
idx = towlib.find('def fold_feed')
print(towlib[idx:idx+3500])
