import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
import mpmath as mp

mp.mp.dps = 150

# Numerical reference
def depth2(x, y):
    """Li_{2,1}(x,y) = sum_{n>m>=1} x^n y^m / (n^2 m)."""
    rho = max(abs(x), abs(x * y))
    if not 0 < rho < 1:
        raise ValueError(f"The series requires |x|, |xy| < 1")
    target = mp.power(10, -mp.mp.dps - 5)
    nmax = int(mp.ceil(mp.log(target * (1 - rho)) / mp.log(rho)))
    h, xy_power, total = mp.mpf(0), x * y, mp.mpf(0)
    for n in range(1, nmax + 1):
        total += h / n**2
        h = x * h + x * xy_power / n
        xy_power *= x * y
    return total

r2 = mp.sqrt(2)
lam_mp = depth2(1 - r2, -r2/2)
swap_mp = depth2(1 - r2/2, -r2)

print("="*70)
print("Building depth-2 H-word vectors")
print("="*70)

print(f"\nLam = Li_{{2,1}}(1-sqrt(2), -sqrt(2)/2):")
print(f"  Numerical value: {mp.nstr(lam_mp, 40)}")

print(f"\nswap = Li_{{2,1}}(1-sqrt(2)/2, -sqrt(2)):")
print(f"  Numerical value: {mp.nstr(swap_mp, 40)}")

# ========== Lam = Li_{2,1}(x1, y1) where x1 = 1-sqrt(2), y1 = -sqrt(2)/2 ==========
# Letters: a1 = 1/x1 = -A, a2 = 1/(x1*y1) 
# x1*y1 = (1-sqrt(2)) * (-sqrt(2)/2) = -sqrt(2)/2 + 1 = 1 - sqrt(2)/2
# 1/(x1*y1) = 1/(1 - sqrt(2)/2) = 2/(2-sqrt(2)) = 2(2+sqrt(2))/2 = 2+sqrt(2)

# So: a1 = -A (in POINTS as -A)
# a2 = 2+sqrt(2) = A+1 (cross ratio from frame)

# G(0, a1, a2; 1) = iterated integral from 0 to 1 with letters a1, a2

# Find a good frame: we need p0, p1, pinf such that:
# 1. The endpoint is 1 (integral from 0 to 1)
# 2. The letters a1 = -A and a2 = A+1 are both in cross-ratio form relative to the frame

# The simplest is to work backwards: we want G(0, a1, a2; 1)
# In forms_path, this corresponds to integrating forms [dlog(0), dlog(a1), dlog(a2)] from some point to another

# Actually, need to understand: G(a0, a1, a2; z) = int_0^z omega_{a0} int_0^t omega_{a1} int_0^s omega_{a2}
# where omega_a = dt/(t-a) = dlog(t-a)

# Wait, let me reconsider:
# Li_{2,1}(x, y) = sum_{n>m>=1} x^n y^m / (n^2 m)
# = sum_n x^n / n^2 * sum_{m<n} y^m / m
# = sum_n x^n / n^2 * (Li_1(y) - sum_{m>n} y^m / m)
# This doesn't simplify nicely

# Better: Li_{2,1}(x,y) = integral representation
# From verify.py comments: G(0, 1/x, 1/(xy); 1) 

# Actually, I think the G integral uses the Goncharov notation:
# G(a1, a2, a3; z) = ∫_0^z dt1/(t-a1) ∫_0^{t1} dt2/(t2-a2) ∫_0^{t2} dt3/(t3-a3)

# Need to verify the sign. Let me check:
# For Li_3(z) = G(0, 0, 1/z; 1) perhaps?

# Let's just try building the forms and comparing numerically
print("\n" + "="*70)
print("Attempting to build G(0, a1, a2; 1) via forms_path")
print("="*70)

# Get frames for a1 = -A and a2 = A+1
a1 = -C.A
a2 = C.Q(2) + C.A

frame_a1 = C.find_cross(a1)
frame_a2 = C.find_cross(a2)

print(f"\na1 = -A: cross-ratio frame = {frame_a1}")
print(f"a2 = 2+sqrt(2): cross-ratio frame = {frame_a2}")

# Having these frames, we can compute the dlog forms
# But to use forms_path, we need a common frame for integration

# Key insight: We need to integrate from 0 to 1 in the original coordinate
# Pull back via the Möbius transformation

print("\n" + "="*70)
print("Alternative: Construct via shuffle and known paths")
print("="*70)

# Actually, depth-2 hyperlogarithms can be expressed as (Shuffle products)
# Li_{2,1}(x,y) = ... after some manipulation

# For now, let's just compute the DIFF using the reduction expressions
# and check if it lies in the known span