base = '<project>/research/08_theory/'
# find TOW2's prove_E3b34
import os
d2 = base+'TOW2_tower_all_depths/'
print(os.path.exists(d2+'prove_E3b34.py'))
a = open(d2+'prove_E3b34.py').read()
b = open(base+'TOW3_tower_depth6/prove_E3b56.py').read()
print(a[:3000])
