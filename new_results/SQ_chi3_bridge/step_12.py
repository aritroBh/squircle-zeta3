import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
import inspect

print("Source of forms_path:")
print(inspect.getsource(C.forms_path)[:3000])