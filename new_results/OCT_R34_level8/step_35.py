import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

# Let me try building Li_{2,1} directly using similar construction
# For a depth-N iterated integral with letters at specific points

# The key insight from atom() is:
# Li_n(z) corresponds to forms_path([dlog(p0,pinf)] * (n-1) + [dlog(p1,pinf)], POINTS[p0], POINTS[e])
# where (p0, p1, pinf, e) is the frame for z

# For depth-2 Li_{2,1}(x, y), we need a 3-fold iterated integral
# G(0, a1, a2; 1) would be forms_path([dlog(letters for 0), dlog(a1), dlog(a2)], ...)

# But wait, let me check the shuffle approach instead
# Since both sides of the reduction are linear combinations of: 
# Li_3, Li_2 * log, log^3, zeta(2) * log, zeta(3)

# We can compute the DIFF vector directly using the reduction expressions
# and check membership

# Let me directly compute the DIFF for Lam and swap

print("="*70)
print("STEP 2 (alternative): Compute DIFF vectors directly from expressions")
print("="*70)

# Load the reduction expressions
import pickle
reductions = pickle.load(open('research/08_theory/OCT6_octagon_closed/reductions.pkl', 'rb'))

lam_expr = reductions['Lam']
swap_expr = reductions['swap']

print(f"\nLam reduction: {lam_expr}")
print(f"\nswap reduction: {swap_expr}")

# Parse the expressions into vectors
# The symbols used: l2, lu, z2, z3, L2x1, L2x2, L2k, L2d, L2r2h, L2mr2h, L3x1, L3x2, L3k, L3d, L3r2h, L3mr2h

# Define all needed atoms
from fractions import Fraction as Fr

# Li vectors (weight-3 for Li_3, weight-2 for Li_2)
l2 = C.atom(('log', C.Q(2)))  # weight-1
lu = C.atom(('log', C.A))      # weight-1
pi2 = C.atom(('pi2',))
z2 = {w: Fr(c, 6) for w, c in pi2.items()}  # weight-2
z3 = C.atom(('z3',))  # weight-3

# Li_2 atoms (weight-2)
L2x1 = C.atom(('li', 2, C.Q(1) - C.A))
L2x2 = C.atom(('li', 2, C.Q(1) - C.A/2))
L2k = C.atom(('li', 2, C.A - C.Q(2)))
L2d = C.atom(('li', 2, C.Q(2) - C.A))
L2r2h = C.atom(('li', 2, C.A / C.Q(2)))

# Li_2(-sqrt(2)/2) via duplication
L2_half = C.atom(('li', 2, C.Q(1)/2))
L2mr2h = C.add((L2_half, Fr(1,2)), (L2r2h, -1))

print("\nAll L2 atoms defined")

# Li_3 atoms (weight-3)
L3x1 = C.atom(('li', 3, C.Q(1) - C.A))
L3x2 = C.atom(('li', 3, C.Q(1) - C.A/2))
L3k = C.atom(('li', 3, C.A - C.Q(2)))
L3d = C.atom(('li', 3, C.Q(2) - C.A))
L3r2h = C.atom(('li', 3, C.A / C.Q(2)))

# For L3mr2h, need to handle separately
mr2h_cmzv = -C.A / C.Q(2)
try:
    L3mr2h = C.atom(('li', 3, mr2h_cmzv))
    print("L3mr2h from cmzv: SUCCESS")
except:
    print("L3mr2h from cmzv: FAILED, will construct alternative")
    L3mr2h = None

print("\nLi_3 atoms defined")