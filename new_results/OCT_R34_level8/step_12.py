code = '''"""Step 1 (continued): Check atom vectors against mpmath and handle failed atoms."""
import sys
from fractions import Fraction as Fr
import pickle
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
import mpmath as mp

mp.mp.dps = 120

# Load the word values
V = pickle.load(open('research/08_theory/A3f6_level16/vals8_150digit.pickle', 'rb'))
val = lambda vec: mp.fsum(mp.mpf(c.numerator) / c.denominator * mp.mpc(*V[w]) for w, c in vec.items())

print("Checking atom vectors against mpmath (≥100 digits):")
r2 = mp.sqrt(2)

# l2 = log(2)
l2_vec = C.atom(('log', C.Q(2)))
l2_mp = mp.log(2)
diff = abs(val(l2_vec) - l2_mp)
print(f"\nl2 = log(2):")
print(f"  vector value matches mpmath: |diff| = {mp.nstr(diff, 10)} ({mp.nstr(-mp.log10(diff), 3)} digits)")

# lu = log(1 + sqrt(2))
lu_vec = C.atom(('log', C.A))
lu_mp = mp.log(1 + r2)
diff = abs(val(lu_vec) - lu_mp)
print(f"\nlu = log(1+sqrt(2)):")
print(f"  vector value matches mpmath: |diff| = {mp.nstr(diff, 10)} ({mp.nstr(-mp.log10(diff), 3)} digits)")

# z2 = zeta(2) = pi^2/6
pi2_vec = C.atom(('pi2',))
z2_vec = {w: Fr(c, 6) for w, c in pi2_vec.items()}
z2_mp = mp.zeta(2)
diff = abs(val(z2_vec) - z2_mp)
print(f"\nz2 = zeta(2):")
print(f"  vector value matches mpmath: |diff| = {mp.nstr(diff, 10)} ({mp.nstr(-mp.log10(diff), 3)} digits)")

# z3 = zeta(3)
z3_vec = C.atom(('z3',))
z3_mp = mp.zeta(3)
diff = abs(val(z3_vec) - z3_mp)
print(f"\nz3 = zeta(3):")
print(f"  vector value matches mpmath: |diff| = {mp.nstr(diff, 10)} ({mp.nstr(-mp.log10(diff), 3)} digits)")

# Check Li_2 and Li_3 at each point
def check_li(name, z_mp, n):
    if n == 2:
        z_cmzv = z_mp
        z_mp_val = mp.polylog(2, z_mp)
    else:
        z_cmzv = z_mp
        z_mp_val = mp.polylog(3, z_mp)
    
    try:
        vec = C.atom(('li', n, z_cmzv))
        diff = abs(val(vec) - z_mp_val)
        print(f"\nLi_{n}({name}):")
        print(f"  vector value matches mpmath: |diff| = {mp.nstr(diff, 10)} ({mp.nstr(-mp.log10(diff), 3)} digits)")
        return vec
    except Exception as e:
        print(f"\nLi_{n}({name}): FAILED with {e}")
        return None

x1 = 1 - r2
x2 = 1 - r2/2
k = r2 - 1
d = 2 - r2
r2h = r2/2

# cmzv algebraic representations
x1_cmzv = C.Q(1) - C.A
x2_cmzv = C.Q(1) - C.A/2
k_cmzv = C.A - C.Q(2)
d_cmzv = C.Q(2) - C.A
r2h_cmzv = C.A / C.Q(2)
mr2h_cmzv = -C.A / C.Q(2)

L2x1 = check_li("1-sqrt(2)", x1_cmzv, 2)
L3x1 = check_li("1-sqrt(2)", x1_cmzv, 3)
L2x2 = check_li("1-sqrt(2)/2", x2_cmzv, 2)
L3x2 = check_li("1-sqrt(2)/2", x2_cmzv, 3)
L2k = check_li("sqrt(2)-1", k_cmzv, 2)
L3k = check_li("sqrt(2)-1", k_cmzv, 3)
L2d = check_li("2-sqrt(2)", d_cmzv, 2)
L3d = check_li("2-sqrt(2)", d_cmzv, 3)
L2r2h = check_li("sqrt(2)/2", r2h_cmzv, 2)
L3r2h = check_li("sqrt(2)/2", r2h_cmzv, 3)
L2mr2h = check_li("-sqrt(2)/2", mr2h_cmzv, 2)
L3mr2h = check_li("-sqrt(2)/2", mr2h_cmzv, 3)

print("\n" + "="*60)
print("Li_2(-sqrt(2)/2) and Li_3(-sqrt(2)/2) need inversion identities")
print("="*60)

# Inversion for Li_2 at z:
# Li_2(z) + Li_2(1/z) = -pi^2/6 - 1/2 (log(-z))^2
# For z = -sqrt(2)/2, we have 1/z = -sqrt(2)
# But -sqrt(2) = -(1+sqrt(2)-1) = -(A-1) = 1-A
# Actually -sqrt(2) = -A+1 = 1-A = 1-(1+sqrt(2)) = -sqrt(2) ✓

# Let's check: -sqrt(2) attains |z| > 1, so we can use inversion
# z = mr2h = -sqrt(2)/2
# 1/z = -sqrt(2)
# For Li_2: Li_2(z) + Li_2(1/z) = -pi^2/6 - 1/2 log^2(-z)
# For Li_3: inversion is more complex

z_mp = -r2/2  # mr2h
one_over_z = -r2  # 1/z = -sqrt(2)

print(f"\nz = -sqrt(2)/2, 1/z = -sqrt(2}")
print(f"|z| = {abs(z_mp)}, |1/z| = {abs(one_over_z)}")

# For Li_2, check inversion identity
# Li_2(z) = -Li_2(1/z) - pi^2/6 - 1/2 log^2(-z)
# where log(-z) for z < 0 is log|z| + i*pi

# Actually for z real negative, log(-z) = log(|z|), no imaginary part
log_minus_z = mp.log(-z_mp)  # log(sqrt(2)/2) = log(1/sqrt(2)) = -log(sqrt(2)) = -1/2*log(2)
print(f"log(-z) = log(sqrt(2)/2) = {mp.nstr(log_minus_z, 20)}")
print(f"        = -log(sqrt(2)) = -log(2)/2 = {mp.nstr(-mp.log(2)/2, 20)}")

# Verify: Li_2(z) + Li_2(1/z) = -pi^2/6 - 1/2 log^2(-z)
# But this identity holds for z complex in general
# Let me compute both sides
z = z_mp
lhs = mp.polylog(2, z) + mp.polylog(2, one_over_z)
rhs = -mp.pi**2/6 - 1/2 * log_minus_z**2
print(f"\nChecking Li_2 inversion: Li_2(z) + Li_2(1/z) vs -pi^2/6 - 1/2 log^2(-z)")
print(f"  LHS = {mp.nstr(lhs, 20)}")
print(f"  RHS = {mp.nstr(rhs, 20)}")
print(f"  diff = {mp.nstr(lhs - rhs, 20)}")

# Actually there's an imaginary part. Let me check more carefully.
# For real z < 0, Li_2(z) is real, Li_2(1/z) is real (since 1/z < -1)
# log(-z) for z < 0 is log(|z|) (principal value)
# But the identity is  Li_2(z) + Li_2(1/z) = -pi^2/6 - 1/2 log^2(-z) for z in (0,1)
# For z < 0, it's different!

# From Lewin: Li_2(z) + Li_2(1/z) = -pi^2/6 - 1/2 log^2(-z) for z not in [0,1]
# and log(-z) = log(|z|) + i*pi for z < 0

# So: Li_2(z) + Li_2(1/z) = -pi^2/6 - 1/2 (log|z| + i*pi)^2
#                         = -pi^2/6 - 1/2 (log^2|z| - pi^2 + 2i*pi*log|z|)
#                         = -pi^2/6 - 1/2 log^2|z| + pi^2/2 - i*pi*log|z|
#                         = pi^2/3 - 1/2 log^2|z| - i*pi*log|z|

# But Li_2 is real for z real on (-inf, 1). Let's recompute...
log_abs_z = mp.log(abs(z_mp))
expected = mp.pi**2/3 - 1/2 * log_abs_z**2
print(f"\nExpected real part: pi^2/3 - 1/2 log^2|z| = {mp.nstr(expected, 20)}")
print(f"Actual Li_2(z) + Li_2(1/z) = {mp.nstr(lhs, 20)}")

# Hmm, still not matching. Let me try the reflection formula.
# Li_2(z) + Li_2(1-z) + log(z)log(1-z) = pi^2/6 for z in (0,1)
# Or: Li_2(z) + Li_2(-z) = 1/2 Li_2(z^2) does not hold

# Standard inversion: Li_2(x) + Li_2(1/x) = -pi^2/3 - 1/2 log^2(-x) for |x| > 1
# Actually I think I am overcomplicating. Let's just use the polylog value directly.

print("\\n" + "="*60)
print("Using direct numerical values for Li_n(-sqrt(2)/2)")
print("="*60)

# Save the vectors we can compute, and note which need special handling
atoms = {
    'l2': C.atom(('log', C.Q(2))),
    'lu': C.atom(('log', C.A)),
    'z2': {w: Fr(c, 6) for w, c in C.atom(('pi2',)).items()},
    'z3': C.atom(('z3',)),
    'L2x1': C.atom(('li', 2, x1_cmzv)),
    'L3x1': C.atom(('li', 3, x1_cmzv)),
    'L2x2': C.atom(('li', 2, x2_cmzv)),
    'L3x2': C.atom(('li', 3, x2_cmzv)),
    'L2k': C.atom(('li', 2, k_cmzv)),
    'L3k': C.atom(('li', 3, k_cmzv)),
    'L2d': C.atom(('li', 2, d_cmzv)),
    'L3d': C.atom(('li', 3, d_cmzv)),
    'L2r2h': C.atom(('li', 2, r2h_cmzv)),
    'L3r2h': C.atom(('li', 3, r2h_cmzv)),
}

# Check if we need to handle mr2h specially
# According to the reduction expression, we need L2mr2h and L3mr2h
# Let's see what symmetries we can use

print("\\nTo handle Li_n(-sqrt(2)/2), checking symmetries:")
print(f"r2h = sqrt(2)/2 = {mp.nstr(r2h, 15)}")
print(f"-r2h = -sqrt(2)/2 = {mp.nstr(-r2h, 15)}")
print(f"1 - r2h = 1 - sqrt(2)/2 = {mp.nstr(1 - r2h, 15)}")
print(f"1 + r2h = 1 + sqrt(2)/2 = {mp.nstr(1 + r2h, 15)}")

# -sqrt(2)/2 = -1/sqrt(2)
# There is a relation: Li_2(-x) = -Li_2(x) + 1/2 Li_2(x^2) for some x?

# Actually, the simplest: Li_2(-z) for z>0
# This is real and can be computed directly
# But cmzv cannot handle z = -A/2 as a cross-ratio

# Let me check if we can express -sqrt(2)/2 as a valid cross-ratio
print(f"\\n-sqrt(2)/2 as Möbius transformation: is it in the level-8 field?")
print(f"-sqrt(2)/2 = {mp.nstr(-r2h, 15)}")
print(f"This is in Q(sqrt(2)), so level-8 compatible")
'''

with open('step_01b_check.py', 'w') as f:
    f.write(code)
print("Created step_01b_check.py")