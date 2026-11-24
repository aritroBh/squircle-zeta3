import os
base = '<project>'
p = os.path.join(base, 'research/08_theory/TOW1_tower_law')
print(os.listdir(p))
with open(os.path.join(p, 'towlib.py')) as f:
    towlib = f.read()
print(towlib)
