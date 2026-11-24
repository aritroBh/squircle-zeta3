# Let me rewrite to use full paths and write locally in work directory
code = '''"""Step 1: Define atoms as H-word vectors and check against mpmath."""
import sys
from fractions import Fraction as Fr
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
import mpmath as mp

mp.mp.dps = 120

# Define all algebraic quantities we need
r2 = mp.sqrt(2)
x1 = 1 - r2  # 1 - sqrt(2)
x2 = 1 - r2/2  # 1 - sqrt(2)/2  
k = r2 - 1  # sqrt(2) - 1
d = 2 - r2  # 2 - sqrt(2)
r2h = r2/2  # sqrt(2)/2
mr2h = -r2/2  # -sqrt(2)/2

print("Algebraic quantities defined:")
print(f"  x1 = 1 - sqrt(2) = {mp.nstr(x1, 20)}")
print(f"  x2 = 1 - sqrt(2)/2 = {mp.nstr(x2, 20)}")
print(f"  k = sqrt(2) - 1 = {mp.nstr(k, 20)}")
print(f"  d = 2 - sqrt(2) = {mp.nstr(d, 20)}")
print(f"  r2h = sqrt(2)/2 = {mp.nstr(r2h, 20)}")
print(f"  mr2h = -sqrt(2)/2 = {mp.nstr(mr2h, 20)}")

# Now define atoms in cmzv algebra
print("\\nDefining atoms as H-word vectors:")

# l2 = log(2)
l2 = C.atom(('log', C.Q(2)))
print(f"l2 = log(2): vector with {len(l2)} words")

# lu = log(1 + sqrt(2))
lu = C.atom(('log', C.A))
print(f"lu = log(1+sqrt(2)): vector with {len(lu)} words")

# z2 = zeta(2) = pi^2/6 - we use ('pi2',) which is pi^2
pi2 = C.atom(('pi2',))
z2_vec = {w: Fr(c, 6) for w, c in pi2.items()}
print(f"z2 = pi^2/6: vector with {len(z2_vec)} words")

# z3 = zeta(3)
z3 = C.atom(('z3',))
print(f"z3 = zeta(3): vector with {len(z3)} words")

# Li_n at algebraic points
# Need to convert each point to Q basis representation

# x1 = 1 - sqrt(2) = 1 - A where A = 1 + sqrt(2)
x1_cmzv = C.Q(1) - C.A
print(f"\\nx1 = 1 - sqrt(2) = {x1_cmzv}")

try:
    L2x1 = C.atom(('li', 2, x1_cmzv))
    print(f"L2x1 = Li_2(1-sqrt(2)): vector with {len(L2x1)} words")
except Exception as e:
    print(f"L2x1 FAILED: {e}")
    L2x1 = None

try:
    L3x1 = C.atom(('li', 3, x1_cmzv))
    print(f"L3x1 = Li_3(1-sqrt(2)): vector with {len(L3x1)} words")
except Exception as e:
    print(f"L3x1 FAILED: {e}")
    L3x1 = None

# x2 = 1 - sqrt(2)/2 = 1 - A/2
x2_cmzv = C.Q(1) - C.A/2
print(f"\\nx2 = 1 - sqrt(2)/2 = {x2_cmzv}")

try:
    L2x2 = C.atom(('li', 2, x2_cmzv))
    print(f"L2x2 = Li_2(1-sqrt(2)/2): vector with {len(L2x2)} words")
except Exception as e:
    print(f"L2x2 FAILED: {e}")
    L2x2 = None

try:
    L3x2 = C.atom(('li', 3, x2_cmzv))
    print(f"L3x2 = Li_3(1-sqrt(2)/2): vector with {len(L3x2)} words")
except Exception as e:
    print(f"L3x2 FAILED: {e}")
    L3x2 = None

# k = sqrt(2) - 1 = A - 2
k_cmzv = C.A - C.Q(2)
print(f"\\nk = sqrt(2) - 1 = {k_cmzv}")

try:
    L2k = C.atom(('li', 2, k_cmzv))
    print(f"L2k = Li_2(sqrt(2)-1): vector with {len(L2k)} words")
except Exception as e:
    print(f"L2k FAILED: {e}")
    L2k = None

try:
    L3k = C.atom(('li', 3, k_cmzv))
    print(f"L3k = Li_3(sqrt(2)-1): vector with {len(L3k)} words")
except Exception as e:
    print(f"L3k FAILED: {e}")
    L3k = None

# d = 2 - sqrt(2) = 2 - A
d_cmzv = C.Q(2) - C.A
print(f"\\nd = 2 - sqrt(2) = {d_cmzv}")

try:
    L2d = C.atom(('li', 2, d_cmzv))
    print(f"L2d = Li_2(2-sqrt(2)): vector with {len(L2d)} words")
except Exception as e:
    print(f"L2d FAILED: {e}")
    L2d = None

try:
    L3d = C.atom(('li', 3, d_cmzv))
    print(f"L3d = Li_3(2-sqrt(2)): vector with {len(L3d)} words")
except Exception as e:
    print(f"L3d FAILED: {e}")
    L3d = None

# r2h = sqrt(2)/2 = A/2
r2h_cmzv = C.A / C.Q(2)
print(f"\\nr2h = sqrt(2)/2 = {r2h_cmzv}")

try:
    L2r2h = C.atom(('li', 2, r2h_cmzv))
    print(f"L2r2h = Li_2(sqrt(2)/2): vector with {len(L2r2h)} words")
except Exception as e:
    print(f"L2r2h FAILED: {e}")
    L2r2h = None

try:
    L3r2h = C.atom(('li', 3, r2h_cmzv))
    print(f"L3r2h = Li_3(sqrt(2)/2): vector with {len(L3r2h)} words")
except Exception as e:
    print(f"L3r2h FAILED: {e}")
    L3r2h = None

# mr2h = -sqrt(2)/2 = -A/2
mr2h_cmzv = -C.A / C.Q(2)
print(f"\\nmr2h = -sqrt(2)/2 = {mr2h_cmzv}")

try:
    L2mr2h = C.atom(('li', 2, mr2h_cmzv))
    print(f"L2mr2h = Li_2(-sqrt(2)/2): vector with {len(L2mr2h)} words")
except Exception as e:
    print(f"L2mr2h FAILED: {e}")
    L2mr2h = None

try:
    L3mr2h = C.atom(('li', 3, mr2h_cmzv))
    print(f"L3mr2h = Li_3(-sqrt(2)/2): vector with {len(L3mr2h)} words")
except Exception as e:
    print(f"L3mr2h FAILED: {e}")
    L3mr2h = None

print("\\n" + "="*60)
print("Step 1 complete: All atoms defined")
'''
with open('step_01_atoms.py', 'w') as f:
    f.write(code)
print("Created step_01_atoms.py in work directory")