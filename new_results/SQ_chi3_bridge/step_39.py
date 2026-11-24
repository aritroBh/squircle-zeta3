import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
from fractions import Fraction
import pickle
import flint
from sympy.ntheory import sqrt_mod

# Load weight3_with_lifts.pickle and k3_newHrows.pickle
w = pickle.load(open("research/08_theory/A3f3_level8_compute/weight3_with_lifts.pickle", "rb"))
k = pickle.load(open("research/08_theory/A3f3_level8_compute/k3_newHrows.pickle", "rb"))

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

# Set up column indexing
cols = list(w['cols'])
idx = {c: i for i, c in enumerate(cols)}
n = len(cols)

def red(x, trunc=False, p=0):
    """Convert a Fraction to a field element mod p."""
    x = Fraction(x)
    if trunc:
        return int(x) % p
    return (x.numerator % p) * pow(x.denominator % p, -1, p) % p

def rank_mod(rows, p):
    """Compute rank of matrix mod p."""
    M = flint.nmod_mat(len(rows), n, [v for r in rows for v in r], p)
    return M.rank()

def vector_to_field(vec, p, s):
    """Convert a vector (dict of word -> coeff) to a list mod p."""
    result = [0] * n
    for word, coeff in vec.items():
        result[idx[word]] = (result[idx[word]] + red(coeff, p=p)) % p
    return result

def krow_to_field(krow, p, s):
    """Convert a krow (word -> (x, y) meaning x + y*√2) to list mod p."""
    result = [0] * n
    for word, (x, y) in krow.items():
        result[idx[word]] = (result[idx[word]] + red(x, p=p) + s * red(y, p=p)) % p
    return result

def check_membership(base_rows, krows, vec, p):
    """Check if vec is in span of (base_rows + krows)."""
    s = sqrt_mod(2, p)
    
    # Build all rows as lists
    W = [vector_to_field(dict(r), p, s) for r in base_rows]
    KR = [krow_to_field(r, p, s) for r in krows]
    V = vector_to_field(vec, p, s)
    
    # Rank of base
    r1 = rank_mod(W + KR, p)
    
    # Rank with vec added
    r2 = rank_mod(W + KR + [V], p)
    
    return r1, r2

# Test at two primes
for p in [1000000033, 1000000241]:
    r1, r2 = check_membership(w['rows'], k['krows'], DIFF, p)
    print(f"p = {p}: rank(base) = {r1}, rank(base + DIFF) = {r2}")
    if r1 == r2:
        print(f"  DIFF is IN the span!")
    else:
        print(f"  DIFF is NOT in the span.")