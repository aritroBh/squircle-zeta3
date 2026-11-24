"""Complete proof script for OCT_R34_level8 reductions R3/R4."""
import sys, pickle, time
from fractions import Fraction as Fr
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
import mpmath as mp
import flint
from sympy.ntheory import sqrt_mod

mp.mp.dps = 150

print("="*70)
print("OCT_R34_level8: Independent exact proof of octagon depth-2 reductions")
print("="*70)

# ====================================================================================
# STEP 1: Define all atoms as H-word vectors
# ====================================================================================
print("\n" + "="*70)
print("STEP 1: Define atoms as H-word vectors")
print("="*70)

# Load word values for numerical checks
V = pickle.load(open('research/08_theory/A3f6_level16/vals8_150digit.pickle', 'rb'))
val = lambda vec: mp.fsum(mp.mpf(c.numerator) / mp.mpf(c.denominator) * mp.mpc(*V[w]) for w, c in vec.items())

r2 = mp.sqrt(2)

# Basic atoms
l2_vec = C.atom(('log', C.Q(2)))  # log(2)
lu_vec = C.atom(('log', C.A))      # log(1 + sqrt(2))
pi2_vec = C.atom(('pi2',))         # pi^2
z2_vec = {w: Fr(c, 6) for w, c in pi2_vec.items()}  # zeta(2) = pi^2/6
z3_vec = C.atom(('z3',))           # zeta(3)

print("Basic atoms: l2, lu, z2, z3 loaded")

# Li_n at algebraic points that cmzv can handle
x1_cmzv = C.Q(1) - C.A  # 1 - sqrt(2)
x2_cmzv = C.Q(1) - C.A/2  # 1 - sqrt(2)/2
k_cmzv = C.A - C.Q(2)   # sqrt(2) - 1
d_cmzv = C.Q(2) - C.A   # 2 - sqrt(2)
r2h_cmzv = C.A / C.Q(2)  # sqrt(2)/2

L2x1 = C.atom(('li', 2, x1_cmzv))
L3x1 = C.atom(('li', 3, x1_cmzv))
L2x2 = C.atom(('li', 2, x2_cmzv))
L3x2 = C.atom(('li', 3, x2_cmzv))
L2k = C.atom(('li', 2, k_cmzv))
L3k = C.atom(('li', 3, k_cmzv))
L2d = C.atom(('li', 2, d_cmzv))
L3d = C.atom(('li', 3, d_cmzv))
L2r2h = C.atom(('li', 2, r2h_cmzv))
L3r2h = C.atom(('li', 3, r2h_cmzv))

print("Li_n vectors at x1, x2, k, d, r2h loaded")

# For L2mr2h, use duplication formula: Li_2(-z) = 1/2 Li_2(z^2) - Li_2(z)
L2_half = C.atom(('li', 2, C.Q(1)/2))  # Li_2(1/2)
L2mr2h = C.add((L2_half, Fr(1,2)), (L2r2h, -1))

print("L2mr2h from duplication: Li_2(-sqrt(2)/2) = 1/2 Li_2(1/2) - Li_2(sqrt(2)/2)")
print(f"  Checking: |vec - mpmath| = {mp.nstr(abs(val(L2mr2h) - mp.polylog(2, -r2/2)), 10)}")

# For L3mr2h, check if direct atom works
mr2h_cmzv = -C.A / C.Q(2)
try:
    L3mr2h = C.atom(('li', 3, mr2h_cmzv))
    print(f"L3mr2h directly from cmzv: SUCCESS")
    print(f"  Checking: |vec - mpmath| = {mp.nstr(abs(val(L3mr2h) - mp.polylog(3, -r2/2)), 10)}")
except ValueError as e:
    print(f"L3mr2h directly: FAILED ({e})")
    L3mr2h = None

# ====================================================================================
# Numerical reference for depth-2 values
# ====================================================================================
def depth2(x, y):
    """Li_{2,1}(x,y) = sum_{n>m>=1} x^n y^m / (n^2 m)."""
    rho = max(abs(x), abs(x * y))
    if not 0 < rho < 1:
        raise ValueError(f"The series requires |x|, |xy| < 1, got rho = {rho}")
    target = mp.power(10, -mp.mp.dps - 5)
    nmax = int(mp.ceil(mp.log(target * (1 - rho)) / mp.log(rho)))
    h, xy_power, total = mp.mpf(0), x * y, mp.mpf(0)
    for n in range(1, nmax + 1):
        total += h / n**2
        h = x * h + x * xy_power / n
        xy_power *= x * y
    return total

lam_mp = depth2(1 - r2, -r2/2)
swap_mp = depth2(1 - r2/2, -r2)

print(f"\nNumerical reference values:")
print(f"  Lam = Li_{{2,1}}(1-sqrt(2), -sqrt(2)/2) = {mp.nstr(lam_mp, 30)}")
print(f"  swap = Li_{{2,1}}(1-sqrt(2)/2, -sqrt(2)) = {mp.nstr(swap_mp, 30)}")

print("\nSTEP 1 complete")
