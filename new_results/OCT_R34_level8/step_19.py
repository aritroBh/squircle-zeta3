Now I understand the structure better. Let me write a complete script for steps 1-5.

code = r'''"""Complete proof script for OCT_R34_level8 reductions R3/R4."""
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
val = lambda vec: mp.fsum(mp.mpf(c.numerator) / c.denominator * mp.mpc(*V[w]) for w, c in vec.items())

r2 = mp.sqrt(2)

# Basic atoms
l2_vec = C.atom(('log', C.Q(2)))  # log(2)
lu_vec = C.atom(('log', C.A))      # log(1 + sqrt(2))
pi2_vec = C.atom(('pi2',))         # pi^2
z2_vec = {w: Fr(c, 6) for w, c in pi2_vec.items()}  # zeta(2) = pi^2/6
z3_vec = C.atom(('z3',))           # zeta(3)

print(f"Basic atoms: l2, lu, z2, z3 loaded")

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

print(f"Li_n vectors at x1, x2, k, d, r2h loaded")

# For L2mr2h and L3mr2h (Li_n at -sqrt(2)/2), use duplication formula
# Li_2(-z) = 1/2 Li_2(z^2) - Li_2(z) for z > 0
# For z = sqrt(2)/2, we have z^2 = 1/2

L2_half = C.atom(('li', 2, C.Q(1)/2))  # Li_2(1/2)
L2mr2h = C.add((L2_half, Fr(1,2)), (L2r2h, -1))  # using duplication

print(f"L2mr2h from duplication: Li_2(-sqrt(2)/2) = 1/2 Li_2(1/2) - Li_2(sqrt(2)/2)")

# Check L2mr2h numerically
L2mr2h_mp = mp.polylog(2, -r2/2)
print(f"  Checking: |vec - mpmath| = {mp.nstr(abs(val(L2mr2h) - L2mr2h_mp), 10)}")

# For L3mr2h, we don't have a simple formula. Check if -A/2 works...
mr2h_cmzv = -C.A / C.Q(2)
try:
    L3mr2h = C.atom(('li', 3, mr2h_cmzv))
    print(f"L3mr2h directly from cmzv: SUCCESS")
except ValueError as e:
    print(f"L3mr2h directly: FAILED ({e})")
    # Try to construct via identities
    # For now, we'll handle this case specially
    print(f"  L3mr2h will be handled through span membership")
    # Create placeholder - this needs proper handling
    L3mr2h = None

# ====================================================================================
# STEP 2: Build depth-2 vectors for Li_{2,1}
# ====================================================================================
print("\n" + "="*70)
print("STEP 2: Depth-2 vectors for Lam and swap")
print("="*70)

# Numerical reference function
def depth2(x, y, p=2, details=False):
    """Sum Li_{p,1} = sum_{n>m>=1} x^n y^m / (n^p m)."""
    if p not in (1, 2):
        raise ValueError("Only Li_{1,1} and Li_{2,1} are used here")
    rho = max(abs(x), abs(x * y))
    if not 0 < rho < 1:
        if rho == 0:
            return mp.mpf(0)
        raise ValueError("The defining series requires |x|, |xy| < 1")
    target = mp.power(10, -mp.mp.dps - 5)
    nmax = int(mp.ceil(mp.log(target * (1 - rho)) / mp.log(rho)))
    h, xy_power, total = mp.mpf(0), x * y, mp.mpf(0)
    for n in range(1, nmax + 1):
        total += h / n**p
        h = x * h + x * xy_power / n
        xy_power *= x * y
    return (total, nmax, rho**(nmax+1)/(1-rho)) if details else total

# Compute numerical values
x1_mp = 1 - r2
x2_mp = 1 - r2/2
lam_mp = depth2(x1_mp, -r2/2)
swap_mp = depth2(x2_mp, -r2)

print(f"\nNumerical reference values:")
print(f"  Lam = Li_{{2,1}}(1-sqrt(2), -sqrt(2)/2) = {mp.nstr(lam_mp, 30)}")
print(f"  swap = Li_{{2,1}}(1-sqrt(2)/2, -sqrt(2)) = {mp.nstr(swap_mp, 30)}")

# Now build depth-2 vectors using forms_path
# The letters for Lam: 1/x = -(1+sqrt(2)) = -A, 1/(xy) = 2+sqrt(2)
# The letters for swap: 1/x = 2+sqrt(2), 1/(xy) = -(1+sqrt(2)) = -A

# Note: C.POINTS = {None:0, Q(0):4, Q(1):6, Q(-1):2, A:7, B:5, -A:1, -B:3}
# We have -A in the configuration!
# And 2+sqrt(2) = A+1 which we need to check

# Let's see if 2+sqrt(2) is accessible
# 2+sqrt(2) = (1+sqrt(2)) + 1 = A + 1
# Cross ratio: cross(e, p0, p1, pinf)

# Try to find if A+1 is a valid cross-ratio
# From the POINTS dict, we have: 0, 1, -1, A, B, -A, -B
# B = 1 - A = -sqrt(2), so -B = sqrt(2)
# A+1 = 2+sqrt(2) is NOT directly in the list

# But we can compute it via cross ratio formula
# cross(e, p0, p1, pinf) = (e-p0)/(e-pinf) * (p1-pinf)/(p1-p0)

# Let me check: for a subset of points, what cross ratios can we get?
print("\nEnumerating cross ratios from configuration points...")
points = [C.Q(0), C.Q(1), C.Q(-1), C.A, C.B, -C.A, -C.B]
cross_ratios = set()
for i, e in enumerate(points):
    for j, p0 in enumerate(points):
        for k, p1 in enumerate(points):
            for l, pinf in enumerate(points):
                if len(set([i,j,k,l])) == 4:
                    cr = C.cross(e, p0, p1, pinf)
                    cross_ratios.add(cr)

# Check which values we need
a1_Lam = C.Q(-1) - C.A  # -A
a2_Lam = C.Q(2) + C.A    # 2 + sqrt(2)

print(f"\nLam letters:")
print(f"  a1 = -A = -(1+sqrt(2)) = {a1_Lam}")
print(f"  a2 = 2+sqrt(2) = {a2_Lam}")
print(f"  -A in cross ratios: {a1_Lam in cross_ratios}")
print(f"  2+sqrt(2) in cross ratios: {a2_Lam in cross_ratios}")

# Find a Möbius frame
# The suggestion is: g takes (0, 1, ∞) to (0, 1, -1-sqrt(2))
# which sends -A to -1 and 2+sqrt(2) to ∞

# Actually, let's try to find a working frame by enumeration
print("\nSearching for Möbius frame for Lam letters...")

# The forms are: [dlog(a1, pinf), dlog(a2, pinf)] plus the initial dlog(0, pinf) for G(0, a1, a2; 1)
# G(0, a1, a2; 1) integrated from 0 to 1

# According to the convention in cmzv.atom for Li_n:
# Li_n(z) is built via C.forms_path([C.dlog(a, q), ...], p0, p1) where z = cross(a, p0, p1, q)

# For depth 2: we need
# G(0, a1, a2; 1) = forms_path([dlog(0, pinf), dlog(g(a1), pinf), dlog(g(a2), pinf)], p0, p1)

# But the G-integral convention may differ. Let me try using H-words directly through path expansion

# Actually, cmzv may have a built-in for iterated integrals. Let me check for a depth2 function

'''
with open('step_03_full.py', 'w') as f:
    f.write(code)
print("Created step_03_full.py (part 1)")