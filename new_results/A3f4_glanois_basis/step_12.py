q = 'research/08_theory/A3f4_glanois_basis/Q/'
import os
print(os.listdir(q))
src = open(q+'cmzv.py').read()
print(len(src))
print(src[:5000])
