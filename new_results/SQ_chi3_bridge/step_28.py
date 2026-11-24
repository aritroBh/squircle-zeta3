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

# Check all values
print("Checking important values against mpmath:\n")

# X = Li3(b)
X = C.atom(('li', 3, C.B))
X_re, X_im = compute_vec_value(X, vals)
b_val = mp.sqrt(2) - 1
X_mp = mp.polylog(3, b_val)
print(f"Li3(b):")
print(f"  Engine: {X_re}")
print(f"  mpmath: {X_mp}")
print(f"  Diff: {abs(X_re - X_mp)}")
print(f"  Digits: {-mp.log(abs(X_re - X_mp))/mp.log(10)}")

# Y = Li3(-b)
Y = C.atom(('li', 3, -C.B))
Y_re, Y_im = compute_vec_value(Y, vals)
Y_mp = mp.polylog(3, -b_val)
print(f"\nLi3(-b):")
print(f"  Engine: {Y_re}")
print(f"  mpmath: {Y_mp}")
print(f"  Diff: {abs(Y_re - Y_mp)}")
print(f"  Digits: {-mp.log(abs(Y_re - Y_mp))/mp.log(10)}")

# Lv = log(1+√2)
Lv = C.atom(('log', C.A))
Lv_re, Lv_im = compute_vec_value(Lv, vals)
Lv_mp = mp.log(1 + mp.sqrt(2))
print(f"\nlog(1+√2):")
print(f"  Engine: {Lv_re}")
print(f"  mpmath: {Lv_mp}")
print(f"  Diff: {abs(Lv_re - Lv_mp)}")
print(f"  Digits: {-mp.log(abs(Lv_re - Lv_mp))/mp.log(10)}")

# P2 = π²
P2 = C.atom(('pi2',))
P2_re, P2_im = compute_vec_value(P2, vals)
P2_mp = mp.pi**2
print(f"\nπ²:")
print(f"  Engine: {P2_re}")
print(f"  mpmath: {P2_mp}")
print(f"  Diff: {abs(P2_re - P2_mp)}")
print(f"  Digits: {-mp.log(abs(P2_re - P2_mp))/mp.log(10)}")