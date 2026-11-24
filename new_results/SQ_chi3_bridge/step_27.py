import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
from fractions import Fraction
import pickle
import mpmath as mp

mp.mp.prec = 500

# First check X against Li3(b)
X = C.atom(('li', 3, C.B))

# Load vals
vals_path = 'research/08_theory/A3f6_level16/vals8_150digit.pickle'
with open(vals_path, 'rb') as f:
    vals = pickle.load(f)

# Compute X value from dictionary
def compute_vec_value(vec, vals):
    total_re = mp.mpf('0')
    total_im = mp.mpf('0')
    for w, coeff in vec.items():
        re_str, im_str = vals[w]
        re_val = mp.mpf(re_str)
        im_val = mp.mpf(im_str)
        total_re += coeff * re_val
        total_im += coeff * im_val
    return total_re, total_im

X_re, X_im = compute_vec_value(X, vals)
print(f"X from engine: {X_re} + {X_im}i")

# Compute Li3(b) from mpmath
b_val = mp.sqrt(2) - 1
X_mpmath = mp.polylog(3, b_val)
print(f"Li3(b) from mpmath: {X_mpmath}")

print(f"\nDifference: {abs(X_re - X_mpmath)}")
print(f"Difference digits: {-mp.log(abs(X_re - X_mpmath))/mp.log(10)}")