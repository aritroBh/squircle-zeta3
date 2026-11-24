import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
import pickle

# Load vals
vals_path = 'research/08_theory/A3f6_level16/vals8_150digit.pickle'
with open(vals_path, 'rb') as f:
    vals = pickle.load(f)

# Check what dlog returns
import inspect
print("Source of dlog:")
print(inspect.getsource(C.dlog))