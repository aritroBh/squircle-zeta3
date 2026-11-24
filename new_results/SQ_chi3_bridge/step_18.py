import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

# Looking at H indices:
# -1 corresponds to dt/t (index -1)
# 0 corresponds to depth regularization

# From POINTS mapping:
# Point 0 -> index 4
# Point 1 -> index 6
# Point -1 -> index 2

# For the zeta values, we need indices corresponding to powers of zeta = e^{iπ/4}
# From the numeric values, we see:
# (-1, -1, 1): has complex value ~0.6635907 + 0.8478278i
# (-1, -1, 7): has complex value ~0.6635907 - 0.8478278i (conjugate)
# These are Li3(zeta) and Li3(zeta^7) = Li3(zeta^{-1})

# The structure is:
# Letter -1: dt/t
# Letter j: ξ^j dt/(1 - ξ^j t)

# For the RHS of the identity:
# 16 * [H(-1, -1, 1) + H(-1, -1, 7) - H(-1, -1, 3) - H(-1, -1, 5)]

# These are depth-1 words after regularization
# Let me check if these exist in the base

# From the task, the words are H(-1, -1, j) for j = 1, 3, 5, 7
# But wait, these should be regularized?

# Let's check what regularize((-1, -1, j)) returns
for j in [0, 1, 2, 3, 4, 5, 6, 7]:
    word = (-1, -1, j)
    reg_word = C.reg(word)
    print(f"reg({word}) = {reg_word}")