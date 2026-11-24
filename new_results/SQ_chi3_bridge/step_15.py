import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
import mpmath as mp

# Now let's set up our variables:
X = C.atom(('li', 3, C.B))
Y = C.atom(('li', 3, -C.B))
Lv = C.atom(('log', C.A))
P2 = C.atom(('pi2',))

# Adjust mpmath precision after importing the engine
mp.mp.prec = 200

# Compute numerical values from mpmath:
b_val = mp.sqrt(2) - 1
L_val = mp.log(1 + mp.sqrt(2))

X_mpmath = mp.polylog(3, b_val)
Y_mpmath = mp.polylog(3, -b_val)
Lv_mpmath = L_val
P2_mpmath = mp.pi**2

print(f"b = √2 - 1 = {b_val}")
print(f"L = log(1+√2) = {L_val}")
print(f"Li3(b) = {X_mpmath}")
print(f"Li3(-b) = {Y_mpmath}")
print(f"π² = {P2_mpmath}")

# Now we need to build the LHS and RHS of the identity:
# LHS = 48 * (X - Y) / 2 = 24 * (X - Y)
# RHS = 16*[H(-1,-1,1) + H(-1,-1,7) - H(-1,-1,3) - H(-1,-1,5)] - 3*P2*Lv + 4*Lv^3

# The H-words are depth-3, with first two indices -1 (dt/t forms)
# Need to construct these as vectors

from functools import reduce

# First, let me understand the structure by looking at regularize and the F vector class
# Let's look at what "reg" function does
import inspect
print("\nSource of reg function:")
print(inspect.getsource(C.reg))