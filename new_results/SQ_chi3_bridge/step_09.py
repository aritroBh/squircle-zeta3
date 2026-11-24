import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
import inspect

# Let's look at the source code for atom to understand the H-word mapping
print("Source of atom function:")
print(inspect.getsource(C.atom)[:2000])