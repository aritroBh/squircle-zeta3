import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
from fractions import Fraction
import pickle
import mpmath as mp

mp.mp.prec = 500

# Load vals
vals_path = 'research/08_theory/A3f6_level16/vals8_150digit.pickle'
with open(vals_path, 'rb') as f:
    vals = pickle.load(f)

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

# Build DIFF
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

DIFF = C.add((LHS, 1), (RHS, -1))

# Compute DIFF value
diff_re, diff_im = compute_vec_value(DIFF, vals)
print(f"DIFF numeric value:")
print(f"  Real: {diff_re}")
print(f"  Imag: {diff_im}")
print(f"  Norm: {mp.sqrt(diff_re**2 + diff_im**2)}")

# Compute digit count
digits = -mp.log(mp.sqrt(diff_re**2 + diff_im**2))/mp.log(10)
print(f"  Digits: {digits}")

# Now let's also check the original identity in mpmath directly
print("\n\nDirect mpmath check of the identity:")
b_val = mp.sqrt(2) - 1
L_val = mp.log(1 + mp.sqrt(2))

# LHS: 48 * chi3(b) = 48 * (Li3(b) - Li3(-b)) / 2 = 24 * (Li3(b) - Li3(-b))
LHS_mp = 24 * (mp.polylog(3, b_val) - mp.polylog(3, -b_val))

# RHS: 32 * sqrt(2) * L(χ8, 3) - 3 * π² * L + 4 * L³
# First compute L(χ8, 3) = (1/√2)(Re Li3(ξ) - Re Li3(ξ³))
xi = mp.e**(mp.j * mp.pi / 4)
L_chi8 = (1/mp.sqrt(2)) * (mp.re(mp.polylog(3, xi)) - mp.re(mp.polylog(3, xi**3)))

RHS_mp = 32 * mp.sqrt(2) * L_chi8 - 3 * mp.pi**2 * L_val + 4 * L_val**3

print(f"\nLHS (48·χ3(b)/2 · 24): {LHS_mp}")
print(f"RHS: {RHS_mp}")
print(f"Difference: {LHS_mp - RHS_mp}")
print(f"Digits: {-mp.log(abs(LHS_mp - RHS_mp))/mp.log(10)}")