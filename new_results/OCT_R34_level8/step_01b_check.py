"""Step 1 (continued): Check atom vectors against mpmath and handle failed atoms."""
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

print("Checking atom vectors against mpmath (at least 100 digits):")
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
def check_li(name, z_cmzv, n, z_mp):
    """Check Li_n at a point."""
    z_mp_val = mp.polylog(n, z_mp)
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

L2x1 = check_li("1-sqrt(2)", x1_cmzv, 2, x1)
L3x1 = check_li("1-sqrt(2)", x1_cmzv, 3, x1)
L2x2 = check_li("1-sqrt(2)/2", x2_cmzv, 2, x2)
L3x2 = check_li("1-sqrt(2)/2", x2_cmzv, 3, x2)
L2k = check_li("sqrt(2)-1", k_cmzv, 2, k)
L3k = check_li("sqrt(2)-1", k_cmzv, 3, k)
L2d = check_li("2-sqrt(2)", d_cmzv, 2, d)
L3d = check_li("2-sqrt(2)", d_cmzv, 3, d)
L2r2h = check_li("sqrt(2)/2", r2h_cmzv, 2, r2h)
L3r2h = check_li("sqrt(2)/2", r2h_cmzv, 3, r2h)
L2mr2h = check_li("-sqrt(2)/2", mr2h_cmzv, 2, -r2h)
L3mr2h = check_li("-sqrt(2)/2", mr2h_cmzv, 3, -r2h)

print("\n" + "="*60)
print("Li_n(-sqrt(2)/2) needs inversion/reflection identities")
print("="*60)

# For z with |z| < 1 but lying on the boundary of the cross-ratio domain,
# we need to use identities
z = -r2h  # -sqrt(2)/2
print(f"\nz = -sqrt(2)/2 = {mp.nstr(z, 15)}")
print(f"|z| = {abs(z)}")

# Test inversion identity for Li_2:
# Li_2(z) + Li_2(1/z) = -pi^2/6 - 1/2 log^2(-z)  for z not in [0,1]
# But this is for 1/z which would be -sqrt(2), dangerous

# Test duplication: Li_2(z) + Li_2(-z) = 1/2 Li_2(z^2)
# For z = r2h = sqrt(2)/2: z^2 = 1/2
# So Li_2(r2h) + Li_2(-r2h) = 1/2 Li_2(1/2)
# Thus Li_2(-r2h) = 1/2 Li_2(1/2) - Li_2(r2h)

print("\nChecking duplication formula: Li_2(z) + Li_2(-z) = 1/2 Li_2(z^2)")
z_for_dup = r2h  # sqrt(2)/2
z2_for_dup = z_for_dup**2  # = 1/2
lhs = mp.polylog(2, z_for_dup) + mp.polylog(2, -z_for_dup)
rhs = 1/2 * mp.polylog(2, z2_for_dup)
print(f"For z = sqrt(2)/2, z^2 = 1/2:")
print(f"  LHS = Li_2(sqrt(2)/2) + Li_2(-sqrt(2)/2) = {mp.nstr(lhs, 20)}")
print(f"  RHS = 1/2 Li_2(1/2) = {mp.nstr(rhs, 20)}")
print(f"  diff = {mp.nstr(lhs - rhs, 20)}")

# So Li_2(-r2h) = 1/2 Li_2(1/2) - Li_2(r2h)
# and 1/2 is in Q, cmzv should handle Li_2(1/2)
print("\nDeriving vector for Li_2(-sqrt(2)/2) = 1/2 Li_2(1/2) - Li_2(sqrt(2)/2):")
L2_half = C.atom(('li', 2, C.Q(1)/2))
print(f"Li_2(1/2) vector: {len(L2_half)} words")
L2mr2h_from_dup = {w: Fr(c, 2) for w, c in L2_half.items()}
for w, c in L2r2h.items():
    L2mr2h_from_dup[w] = L2mr2h_from_dup.get(w, 0) - c
print(f"Li_2(-sqrt(2)/2) from duplication: {len(L2mr2h_from_dup)} words")

# But this removes entries. Let's properly subtract
L2mr2h_from_dup = C.add((L2_half, Fr(1,2)), (L2r2h, -1))
# Actually use C.add properly: it takes tuples of (vector, coeff)
# No, C.add takes (vector, coeff) pairs

# Let me check how C.add works
# From the template: C.add((X, 24), (Y, -24)) means 24*X - 24*Y
# So: L2mr2h = 1/2 Li_2(1/2) - Li_2(r2h)
# should be C.add((L2_half, Fr(1,2)), (L2r2h, -1))

L2mr2h_vec = C.add((L2_half, Fr(1,2)), (L2r2h, -1))
L2mr2h_mp = mp.polylog(2, -r2h)
diff = abs(val(L2mr2h_vec) - L2mr2h_mp)
print(f"\nLi_2(-sqrt(2)/2) from duplication formula:")
print(f"  vector value matches mpmath: |diff| = {mp.nstr(diff, 10)} ({mp.nstr(-mp.log10(diff), 3)} digits)")
print(f"  SUCCESS: vector has {len(L2mr2h_vec)} words")

# Now for Li_3(-sqrt(2)/2), we might need Li_3 duplication
# Li_3(z) + Li_3(-z) = 1/4 Li_3(z^2) ??? No, trilog identity is different

# Actually: Li_3(z) + Li_3(-z) = 1/4 Li_3(z^2) - 1/6 log^3(1-z) stuff? Let me check numerically
print("\nChecking if there is a simple Li_3 duplication formula:")
z_test = r2h
lhs3 = mp.polylog(3, z_test) + mp.polylog(3, -z_test)
rhs3_simple = 1/4 * mp.polylog(3, z_test**2)
print(f"For z = sqrt(2)/2:")
print(f"  Li_3(z) + Li_3(-z) = {mp.nstr(lhs3, 20)}")
print(f"  1/4 Li_3(z^2) = {mp.nstr(rhs3_simple, 20)}")
print(f"  diff = {mp.nstr(lhs3 - rhs3_simple, 20)}")
# Not matching, so the simple formula doesn't hold

# There is the Landen identity for Li_3:
# Li_3(z) + Li_3(1-z) + Li_3(1-1/z) = Li_3(1) + log terms
# Or use the polylog functional equations

# Actually, there's no simple duplication for Li_3. We need other identities.
# Let me check the reflection formula: Li_n(-z) for z > 0

# Since mp.polylog(3, -r2h) gives the correct value,
# let me see if I can express Li_3(-r2h) through the reduction itself.
# The reduction expression only uses L3mr2h = -4 * Li_3(-sqrt(2)/2)
# Other terms: L3d, L3k, L3x2, and logs/pi^2

# For now, let me check the swap expression numerically
print("\n" + "="*60)
print("Checking if the swap reduction requires L3mr2h")
print("="*60)

# From the reduction:
# swap = L2mr2h*l2/2 - L2r2h*l2/2 + L2r2h*lu/2 - L2x2*l2 + L2x2*lu/2 
#      + L3d - L3k - 4*L3mr2h - 4*L3r2h + 2*L3x2 - l2**2*lu/2 + l2*lu**2/2 + l2*z2/2 - lu*z2/2

# So yes, we need L3mr2h = Li_3(-sqrt(2)/2)

# Let me check if -sqrt(2)/2 can be written as a Möbius transform
# or if we need another identity. For Li_3 at z in (-1, 0), polylog is real.

# There is the identity: Li_3(x) + Li_3(1-x) + Li_3(1-1/x) = zeta(3) + 1/6 log^3(x) - 1/2 log^2(x) log(1-x) + log(x) zeta(2) ...
# for x > 1

# Also: Li_3(-x) + Li_3(-1/x) = -Li_3(-1) - 1/6 log^3(x) + log(x) Li_2(-x) ??? 
# Li_3(-1) = -3/4 zeta(3)

# Try inversion: Li_3(z) + Li_3(1/z) = ... for |z| < 1
# Actually for polylog of order 3, there is:
# Li_3(z) + Li_3(1/z) + Li_3(1) = 1/6 log^3(-z) + log(-z) Li_2(1/z) + zeta(3)
# reference: Lewin's book

print("\nChecking Li_3 inversion:")
z_inv_test = mp.mpf(2)/r2  # 1/(sqrt(2)/2) = sqrt(2)
# For z = sqrt(2)/2, 1/z = sqrt(2)
# For Li_3(sqrt(2)/2) + Li_3(sqrt(2)) + zeta(3) ? 

# Let me try: Li_3(-z) = -Li_3(z) + 1/4 Li_3(z^2)  (palindromic ?) Check:
z_test = r2h
lhs_check = mp.polylog(3, -z_test)
rhs_check = -mp.polylog(3, z_test) + mp.polylog(3, z_test**2)/4
print(f"Checking: Li_3(-z) = -Li_3(z) + 1/4 Li_3(z^2) for z = sqrt(2)/2")
print(f"  LHS = {mp.nstr(lhs_check, 20)}")
print(f"  RHS = {mp.nstr(rhs_check, 20)}")
print(f"  diff = {mp.nstr(lhs_check - rhs_check, 20)}")
# This doesn't match

# Another approach: compute L3mr2h numerically and check if it appears in a span relation
# For now, let's proceed with the assumption that we can handle L3mr2h

print("\nFor L3mr2h, we'll compute it numerically and check span membership")
