"""Step 2: Build depth-2 vectors for Li_{2,1}."""
import sys
from fractions import Fraction as Fr
import pickle
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
import mpmath as mp

mp.mp.dps = 150

# Load the numerical reference from verify.py
def depth2(x, y, p=2, details=False):
    """Sum Li_{p,1}; tail <= rho**(N+1)/(1-rho), rho=max(|x|,|xy|)."""
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

# Values for Lam and swap
r2 = mp.sqrt(2)
x1 = 1 - r2  # for Lam: y = -r2/2
x2 = 1 - r2/2  # for swap: y = -r2

# Lam = Li_{2,1}(1-sqrt(2), -sqrt(2)/2)
lam_mp = depth2(x1, -r2/2)
print(f"Lam = Li_{{2,1}}(1-sqrt(2), -sqrt(2)/2) = {mp.nstr(lam_mp, 30)}")

# swap = Li_{2,1}(1-sqrt(2)/2, -sqrt(2))
swap_mp = depth2(x2, -r2)
print(f"swap = Li_{{2,1}}(1-sqrt(2)/2, -sqrt(2)) = {mp.nstr(swap_mp, 30)}")

# The formula: Li_{2,1}(x, y) = G(0, 1/x, 1/(xy); 1) up to a sign
# Let's verify the sign convention by computing both sides

# For Lam: x = x1 = 1 - sqrt(2), y = -r2/2
# 1/x = 1/(1 - sqrt(2)) = (1 + sqrt(2))/(1-2) = -(1+sqrt(2)) = -A
# 1/(xy) = 1/((1 - sqrt(2)) * (-sqrt(2)/2)) = 1/((1 - sqrt(2)) * (-sqrt(2)/2))
#        = -2/(sqrt(2)(1 - sqrt(2))) = 2/(sqrt(2)(sqrt(2) - 1))
#        = 2/(2 - sqrt(2)) = 2/(d) where d = 2 - sqrt(2)

# Let's compute more carefully
x_Lam = x1  # 1 - sqrt(2)
y_Lam = -r2/2  # -sqrt(2)/2

a1 = 1/x_Lam
a2 = 1/(x_Lam * y_Lam)

print(f"\nLam: x = {mp.nstr(x_Lam, 15)}, y = {mp.nstr(y_Lam, 15)}")
print(f"  1/x = {mp.nstr(a1, 15)}")
print(f"  1/(xy) = {mp.nstr(a2, 15)}")
print(f"  -A = -{mp.nstr(1 + r2, 15)}")

# Check: 1/x should be -(1 + sqrt(2)) = -A
print(f"  1/x == -(1+sqrt(2))? {abs(a1 - (-1 - r2)) < 1e-50}")

# The G integral: G(0, a1, a2; 1) with letters a1, a2
# Since a1 = -A, we need to check what vectors we can form

print("\n" + "="*60)
print("Examining cmzv forms_path for depth-2 iterated integrals")
print("="*60)

# Let's look at cmzv to understand forms_path
# The key is: cmzv.forms_path(forms, p0, p1)
# where forms is a list of dlog forms, and p0, p1 are cross-ratio points

# Check what C.POINTS contains
print(f"\nC.POINTS keys: {list(C.POINTS.keys())}")
# Note: cmzv uses cross-ratios as keys

# The forms are C.dlog(a, b) where a, b are in the cross-ratio set
# For a G-integral, we need to express the letters as dlog forms

# First, let's check which cross-ratios are valid
print("\nCross-ratios in configuration (cross-ratio field Q(sqrt(2))):")
# The level-8 configuration has points that form cross-ratios
# Let's find which pairs give the required letters

# Letters for Lam: a1 = -A, a2 = 2 + sqrt(2)
# We need: -A and 2 + sqrt(2) to be in the cross-ratio set

print(f"\nLam letters: a1 = -A = -(1+sqrt(2)), a2 = 2 + sqrt(2)")
print(f"  -A = {mp.nstr(-1 - r2, 20)}")
print(f"  2 + sqrt(2) = {mp.nstr(2 + r2, 20)}")

# For the Möbius frame g taking (0,1,∞) to configuration points
# We need g(a1) and g(a2) also in the configuration, and ∞ not between p0 and p1

# The suggested frame: (p0, p1, pinf) = (0, 1, -1-sqrt(2))
# This sends:
#   0 -> 0
#   1 -> 1
#   ∞ -> -1 - sqrt(2)
# The frame g(z) = z*(1 - (-1-sqrt(2))) / (z - (-1-sqrt(2))/(...)) ?
# Actually need to construct properly...

# Simpler: check what cross-ratios cmzv recognizes
# From the template, we know C.Q(n) and C.A are valid

# Let's enumerate frames as cmzv.find_cross does
# First, understand find_cross

print("\n" + "="*60)
print("Understanding find_cross mechanism")
print("="*60)

# Let me examine what find_cross does for a sample input
# It's in A3f4_glanois_basis/Q/cmzv.py
