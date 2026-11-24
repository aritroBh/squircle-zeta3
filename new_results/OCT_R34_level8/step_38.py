Let me reconsider. Looking at the numerical value: DIFF ≈ -0.23 vs Lam ≈ -0.028. These don't match.

Actually, DIFF should be: vector(Lam) - vector(reduction['Lam'])
where reduction['Lam'] is the expression.
So I need to vector-ize the reduction expression and then compute:

DIFF = V(Lam) - V(reduction)

But I don't have V(Lam) yet! I need to build the depth-2 vector for the actual Lam value.

Let me go back to understanding how to build the depth-2 vector properly.

code = r'''"""Build depth-2 Li_{2,1} vectors via iterated integrals."""
import sys, pickle
from fractions import Fraction as Fr
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
import mpmath as mp

mp.mp.dps = 150

print("="*70)
print("Building depth-2 Li_{2,1} vectors")
print("="*70)

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

print(f"\nNumerical reference:")
print(f"Lam = Li_{{2,1}}(1-sqrt(2), -sqrt(2)/2):")
print(f"     = {mp.nstr(lam_mp, 50)}")

# The key identity:
# Li_{2,1}(x, y) = G(0, 1/x, 1/(xy); 1)  (Goncharov notation)
# where G is the iterated integral: G(a1, a2, a3; z) = ∫_0^z dt1/(t-a1) ∫_0^{t1} dt2/(t2-a2) ∫_0^{t2} dt3/(t3-a3)

# Verify numerically: compute G(0, a1, a2; 1) and check vs depth2
# For Lam: x = 1 - sqrt(2), y = -sqrt(2)/2
# a1 = 1/x = -(1+sqrt(2)) = -A  
# a2 = 1/(xy) 
# xy = (1-sqrt(2))*(-sqrt(2)/2) = sqrt(2)/2 * (sqrt(2)-1) = (2 - sqrt(2))/2 = (2-sqrt(2))/2
# 1/(xy) = 2/(2-sqrt(2)) = 2(2+sqrt(2))/2 = 2+sqrt(2) = A+1

x_mp = 1 - r2
y_mp = -r2/2
a1_mp = 1/x_mp  
a2_mp = 1/(x_mp * y_mp)

print(f"\nFor Lam, verify numerically:")
print(f"  x = {mp.nstr(x_mp, 20)}")
print(f"  y = {mp.nstr(y_mp, 20)}")
print(f"  a1 = 1/x = {mp.nstr(a1_mp, 20)}")
print(f"  a2 = 1/(xy) = {mp.nstr(a2_mp, 20)}")
print(f"  -A = {mp.nstr(-(1+r2), 20)}")
print(f"  A+1 = {mp.nstr(2+r2, 20)}")

# Compute G(0, a1, a2; 1) numerically using mpmath
# G(0, a1, a2; 1) = ∫_0^1 dt1/(t1-0) ∫_0^{t1} dt2/(t2-a1) ∫_0^{t2} dt3/(t3-a2)
#                = ∫_0^1 dt1/t1 ∫_0^{t1} dt2/(t2-a1) ∫_0^{t2} dt3/(t3-a2)

# But this integral diverges at t1=0 (log divergence). We need to check the sign/regularization.

# Alternative: direct series expansion
# Li_{2,1}(x, y) = -G(0, 1/x, 1/(xy); 1) + some boundary terms
# The minus sign comes from the change of variables

# Let me try computing G(a1, a2, a3; 1) numerically:
# G(a1, a2, a3; 1) = ∫_0^1 dt1/(t-a1) ∫_0^{t1} dt2/(t2-a2) ∫_0^{t2} dt3/(t3-a3)

# For 0, a1, a2:
def G_0_a1_a2(a1, a2):
    """Compute G(0, a1, a2; 1) numerically."""
    # G(0, a1, a2; 1) = ∫_0^1 dt1/t1 ∫_0^{t1} dt2/(t2-a1) ∫_0^{t2} dt3/(t3-a2)
    
    # Inner-most: ∫_0^{t2} dt3/(t3-a2) = [log(t3-a2)]_0^{t2} = log(t2-a2) - log(-a2)
    # Middle: ∫_0^{t1} dt2/(t2-a1) * [log(t2-a2) - log(-a2)]
    # Outer: ∫_0^1 dt1/t1 * [result]
    
    # This has div at t1=0. Let's regularize.
    # G(0, a1, a2; 1) - G(0, a1, a2; epsilon) ~ log(epsilon) * I(a1, a2)
    # where I(a1, a2) = ∫_0^1 dt2/(t2-a1) * [log(t2-a2) - log(-a2)]  evaluated from epsilon to 1
    
    # Better: use the differential equation approach
    # d/dz G(a1, a2, a3; z) = 1/(z-a1) * G(a2, a3; z)
    # where G(a1, a2; z) = ∫_0^z dt1/(t1-a1) ∫_0^{t1} dt2/(t2-a2)
    
    # For the case G(0, a1, a2; 1):
    # d/dz G(0, a1, a2; z) = 1/z * G(a1, a2; z)
    # G(a1, a2; z) = ∫_0^z dt1/(t1-a1) ∫_0^{t1} dt2/(t2-a2)
    #              = Li_1(1/(a1/a2))? No...
    
    # Actually, let's compute differently:
    # G(a1, a2, a3; 1) using the shuffle relations and known values
    
    # For a1=0: G(0, a2, a3; 1) 
    # d/dz G(0, a2, a3; z) = 1/z * G(a2, a3; z)
    # G(a2, a3; z) = ∫_0^z dt1/(t1-a2) ∫_0^{t1} dt2/(t2-a3)
    #              = another iterated integral we can compute
    
    # For now, let's use direct numerical integration avoiding the singularity
    def inner(t2):
        # ∫_0^{t2} dt3/(t3-a2) = log((t2-a2)/(-a2)) for real integration
        return mp.log((t2 - a2)/(-a2))
    
    def middle(t1):
        # ∫_0^{t1} dt2/(t2-a1) * inner(t2)
        # Use numerical integration
        if t1 > 1e-100:
            return mp.quad(lambda t2: mp.reciprocal(t2 - a1) * inner(t2), [0, t1])
        else:
            return 0
    
    # Outer integral: regularize by subtracting the divergence
    # G(0, a1, a2; 1) = ∫_0^1 dt1/t1 * middle(t1)
    # Near t1=0, middle(t1) ~ const + O(t1), so integrable
    
    result = mp.quad(lambda t1: mp.reciprocal(t1) * middle(t1), [1e-100, 1])
    return result

# Actually the triple integral is slow. Let me check the relation directly.
# The correct formula involves a shuffle relation.

# Shuffle: Li_1(x) * Li_2(y) = Li_3(xy) + Li_{2,1}(x, y) + Li_{2,1}(y, x)
# where Li_1(x) = -log(1-x)

# For our case: 
# -log(1-x1) * Li_2(y1) = Li_3(x1*y1) + Li_{2,1}(x1, y1) + Li_{2,1}(y1, x1)

# Let x1 = 1 - sqrt(2), y1 = -sqrt(2)/2
x1 = x_mp
y1 = y_mp

print(f"\nUsing shuffle: Li_1(x) * Li_2(y) = Li_3(xy) + Li_{{2,1}}(x,y) + Li_{{2,1}}(y,x)")
print(f"x = {mp.nstr(x1, 20)}, y = {mp.nstr(y1, 20)}")

Li1_x = -mp.log(1 - x1)  # Li_1(x) = -log(1-x)
Li2_y = mp.polylog(2, y1)
Li3_xy = mp.polylog(3, x1 * y1)

print(f"\nLi_1(x) = -log(1-x) = {mp.nstr(Li1_x, 30)}")
print(f"Li_2(y) = {mp.nstr(Li2_y, 30)}")
print(f"Li_3(xy) = {mp.nstr(Li3_xy, 30)}")
print(f"Li_1(x) * Li_2(y) = {mp.nstr(Li1_x * Li2_y, 30)}")

# So: Li_{2,1}(x,y) + Li_{2,1}(y,x) = Li_1(x)*Li_2(y) - Li_3(xy)

rhs_shuffle = Li1_x * Li2_y - Li3_xy
print(f"\nLi_{{2,1}}(x,y) + Li_{{2,1}}(y,x) = {mp.nstr(rhs_shuffle, 30)}")

# We need another relation to separate these two terms.
# One option: use the symmetry/antisymmetry properties

# Let's compute Li_{2,1}(y1, x1) directly:
lam2_mp = depth2(y1, x1)
print(f"Li_{{2,1}}(y, x) = {mp.nstr(lam2_mp, 30)}")
print(f"Li_{{2,1}}(x, y) = {mp.nstr(lam_mp, 30)}")
print(f"Sum = {mp.nstr(lam_mp + lam2_mp, 30)}")

diff_check = (lam_mp + lam2_mp) - rhs_shuffle
print(f"Diff from shuffle: {mp.nstr(diff_check, 10)}")

# This verifies the shuffle relation!

# Now we can use this to express Li_{2,1} in terms of products
# But for the H-word vector, we need to understand the G-integral representation

print("\n" + "="*70)
print("The shuffle relation is verified. Now building G-vector.")
print("="*70)

# The G integral: cmzv can compute G(a1, ..., an; z) via forms_path
# The key is that forms_path takes dlog forms and returns the iterated integral

# For G(a1, a2, a3; z) with a1, a2, a3 in the cross-ratio set:
# G(a1, a2, a3; z) = forms_path([dlog_a1, dlog_a2, dlog_a3], start, end)
# where dlog_a = dlog(point_for_a, point_for_inf)

# Let's compute G(-A, A+1, 0; 1) and see if it relates to Lam
# For this, we need to express the integration in a Möbius frame

# The frame p0, p1, pinf determines the start, end, and infinity
# The dlog points need to be cross ratios in this frame

# Coming back to cmzv-specific construction:
# The atom function uses: forms_path([dlog(p0, pinf)] * (n-1) + [dlog(p1, pinf)], POINTS[p0], POINTS[e])
# So for Li_n(z) with z = cross(e, p0, p1, pinf), the first n-1 letters are at p0, last at p1.

# For G(0, a1, a2; 1), we need to find the correct frame and dlog forms.

# Let me try: the integral from 0 to 1 with letters a1, a2
# means: p0 = 0, p1 = 1 (start and end of integration)
# The letter "0" corresponds to the lower endpoint

# But dlog(p0, pinf) represents dt/(t - p0) pulled back through the frame 
# Actually wait - I think I'm overcomplicating.

# Let me just check what the H-word vector for G looks like:
# Pick a1, a2 that are both valid cross ratios and compute G(a1, a2, 0; 1)

# Actually, let's approach differently: 
# Since we've verified the shuffle relation, we can express Li_{2,1}(x,y) as:
# Li_{2,1}(x,y) = Li_1(x) * Li_2(y) - Li_3(xy) - Li_{2,1}(y,x)

# In terms of H-word vectors:
# V(Li_{2,1})(x,y) = V(Li_1(x1)) * V(Li_2(y1)) - V(Li_3(x1*y1)) - V(Li_{2,1}(y1, x1))
# But this still involves Li_{2,1}!

# The key is that one of the Li_{2,1} terms might have a better representation.
# Or we use more relations.

# For the proof, let's check if the reduction DIFF is in the span
# by treating the depth-2 value as a placeholder and checking membership.

# Actually, the problem states that we should construct the G-vector.
# Let me study the path function more carefully.

print("\nAttempting to use forms_path for G integral...")

# Get POINTS for our letters
print(f"\nC.POINTS = {C.POINTS}")

# For Lam: letters a1 = -A, a2 = A+1, and the "0" letter
# -A = -A is in POINTS as exponent 1
# A+1 = A+1, need to check if it is in POINTS

# Actually, C.POINTS are the configuration points, not the letters.
# The dlog form dlog(p, q) pulls back dt/(t-p)/(t-q) to a specific frame.

# I need to trace through the implementation more carefully.
# Let me check how atom builds Li_3:

z_test = C.Q(1) - C.A
p0, p1, pinf, e = C.find_cross(z_test)
print(f"\nFor z = {z_test}:")
print(f"  Frame: p0={p0}, p1={p1}, pinf={pinf}, e={e}")
print(f"  This means z = cross(e, p0, p1, pinf)")

# The Li_3 vector is built as:
# forms_path([dlog(p0, pinf), dlog(p0, pinf), dlog(p1, pinf)], POINTS[p0], POINTS[e])
forms = [C.dlog(p0, pinf), C.dlog(p0, pinf), C.dlog(p1, pinf)]
vec_Li3_test = C.forms_path(forms, C.POINTS[p0], C.POINTS[e])
print(f"  Vector has {len(vec_Li3_test)} words (should be 64 for Li_3)")

# Now for G(0, a1, a2; 1), we want the integral from 0 to 1
# Letter 0: singular at t=0, maybe dlog(p0, pinf) where p0 = 0?
# Letter a1 = -A: dlog(-A, pinf)
# Letter a2 = A+1: dlog(A+1, pinf)

# But we need pinf such that all letters are valid cross ratios and the endpoint is 1.

# The key is: the integral G(a1, a2, a3; z) can be pulled back to the frame
# and computed using forms_path.

# Let me try with a specific frame:
# p0 = 0, p1 = 1, pinf = something, e = something
# We want the endpoint to be 1 (in the original coordinate)
# Pull back: the forms are [dlog(g(a1)), dlog(g(a2)), dlog(g(0))] where g is the Möbius transform

# This is getting complex. Let me check if there's documentation or examples in the codebase.
'''

with open('step_05_depth2_vec.py', 'w') as f:
    f.write(code)
print("Created step_05_depth2_vec.py")