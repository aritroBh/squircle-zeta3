import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
from fractions import Fraction
import pickle
import mpmath as mp

# Rebuild everything
X = C.atom(('li', 3, C.B))
Y = C.atom(('li', 3, -C.B))
X_minus_Y = C.add((X, 1), (Y, -1))
LHS = {k: v * 24 for k, v in X_minus_Y.items()}

H1 = {(-1, -1, 1): Fraction(1)}
H3 = {(-1, -1, 3): Fraction(1)}
H5 = {(-1, -1, 5): Fraction(1)}
H7 = {(-1, -1, 7): Fraction(1)}
H_sum = C.add((H1, 1), (H7, 1), (H3, -1), (H5, -1))
H_part = {k: v * 16 for k, v in H_sum.items()}

P2 = C.atom(('pi2',))
Lv = C.atom(('log', C.A))
P2_Lv = C.mul(P2, Lv)
Lv2 = C.mul(Lv, Lv)
Lv3 = C.mul(Lv2, Lv)
RHS = C.add((H_part, 1), (P2_Lv, -3), (Lv3, 4))

# DIFF = LHS - RHS
DIFF = C.add((LHS, 1), (RHS, -1))

# Load vals
vals_path = 'research/08_theory/A3f6_level16/vals8_150digit.pickle'
with open(vals_path, 'rb') as f:
    vals = pickle.load(f)

# Compute DIFF numeric value
def compute_vec_value(vec, vals):
    """Compute the value of a vector using the values in vals.
    Returns (real_part, imag_part) as mpmath numbers."""
    total_re = mp.mpf('0')
    total_im = mp.mpf('0')
    for w, coeff in vec.items():
        re_str, im_str = vals[w]
        re_val = mp.mpf(re_str)
        im_val = mp.mpf(im_str)
        total_re += coeff * re_val
        total_im += coeff * im_val
    return total_re, total_im

# Set precision
mp.mp.prec = 500

# Compute DIFF value
diff_re, diff_im = compute_vec_value(DIFF, vals)
print(f"DIFF value (from vals):")
print(f"  real: {diff_re}")
print(f"  imag: {diff_im}")

# Check if DIFF is close to zero
print(f"\nNorm of DIFF: {mp.sqrt(diff_re**2 + diff_im**2)}")