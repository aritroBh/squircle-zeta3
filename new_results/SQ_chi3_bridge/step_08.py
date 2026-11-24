print("Example H-word (should be for Li3(zeta) where zeta is a point):")
# Let's see what H(-1, -1, j) looks like for different j
# From the task, j values are 1, 3, 5, 7 (corresponding to powers of zeta)
# Let me try constructing these

# Looking at the POINTS dict, we have numeric indices:
# 0: 4 entries
# 1: 6 entries (point 1)
# -1: 2 entries (point -1)
# etc.

# From the task, the xi_j form of index j in H(-1, -1, j):
# Need to understand how the engine maps these

# Let's see what atom(('li', 3, z)) returns for various z
import cmzv as C
qs2 = C.Q

# Let's look at the source code for atom to understand the H-word mapping
import inspect
print("Source of atom function:")
print(inspect.getsource(C.atom))