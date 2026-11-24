# Generate the complete result.md
import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
from fractions import Fraction
import pickle
import mpmath as mp
import flint
from sympy.ntheory import sqrt_mod

mp.mp.prec = 500

# Build all vectors
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

# Load vals and compute DIFF numeric value
vals = pickle.load(open("research/08_theory/A3f6_level16/vals8_150digit.pickle", "rb"))

def compute_vec_value(vec, vals):
    total_re = mp.mpf('0')
    total_im = mp.mpf('0')
    for w, coeff in vec.items():
        re_str, im_str = vals[w]
        total_re += coeff * mp.mpf(re_str)
        total_im += coeff * mp.mpf(im_str)
    return total_re, total_im

diff_re, diff_im = compute_vec_value(DIFF, vals)

# Compute direct mpmath check
b_val = mp.sqrt(2) - 1
L_val = mp.log(1 + mp.sqrt(2))
xi = mp.e**(mp.j * mp.pi / 4)
L_chi8 = (1/mp.sqrt(2)) * (mp.re(mp.polylog(3, xi)) - mp.re(mp.polylog(3, xi**3)))

LHS_mp = 24 * (mp.polylog(3, b_val) - mp.polylog(3, -b_val))
RHS_mp = 32 * mp.sqrt(2) * L_chi8 - 3 * mp.pi**2 * L_val + 4 * L_val**3
direct_diff = abs(LHS_mp - RHS_mp)

# Rank checks
w = pickle.load(open("research/08_theory/A3f3_level8_compute/weight3_with_lifts.pickle", "rb"))
k = pickle.load(open("research/08_theory/A3f3_level8_compute/k3_newHrows.pickle", "rb"))

cols = list(w['cols'])
idx = {c: i for i, c in enumerate(cols)}
n = len(cols)

def red(x, trunc=False, p=0):
    x = Fraction(x)
    if trunc:
        return int(x) % p
    return (x.numerator % p) * pow(x.denominator % p, -1, p) % p

def rank_mod(rows, p):
    M = flint.nmod_mat(len(rows), n, [v for r in rows for v in r], p)
    return M.rank()

def vector_to_field(vec, p, s):
    result = [0] * n
    for word, coeff in vec.items():
        result[idx[word]] = (result[idx[word]] + red(coeff, p=p)) % p
    return result

def krow_to_field(krow, p, s):
    result = [0] * n
    for word, (x, y) in krow.items():
        result[idx[word]] = (result[idx[word]] + red(x, p=p) + s * red(y, p=p)) % p
    return result

rank_results = {}
for p in [1000000033, 1000000241]:
    s = sqrt_mod(2, p)
    W = [vector_to_field(dict(r), p, s) for r in w['rows']]
    KR = [krow_to_field(r, p, s) for r in k['krows']]
    V = vector_to_field(DIFF, p, s)
    r1 = rank_mod(W + KR, p)
    r2 = rank_mod(W + KR + [V], p)
    rank_results[p] = (r1, r2)

print(f"DIFF words: {len(DIFF)}")
print(f"DIFF norm: {mp.sqrt(diff_re**2 + diff_im**2)}")
print(f"DIFF digits: {-mp.log(mp.sqrt(diff_re**2 + diff_im**2))/mp.log(10)}")
print(f"Direct check digits: {-mp.log(direct_diff)/mp.log(10)}")
print(f"Rank results: {rank_results}")