import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

# Let me understand the depth-1 words for Li3(zeta^j)
# zeta = e^{i pi/4}, so zeta^j = e^{i j pi/4}
# These correspond to points on the unit circle, which from the information are in the Kummer extension

# But wait, looking at the POINTS, we only see 8 points total, which include points in Q(sqrt2)
# The instructions say: "Depth-1 words for Li3(ξ^j): determine whether H(−1,−1,j) is Li3(ξ^j) or Li3(ξ^{−j}) NUMERICALLY"

# So the words for Li3(zeta^j) should be of form (-1, -1, j) where j = 1, 3, 5, 7

# Let's construct these words and see what they look like
# From the structure, (-1, -1, j) means: dt/t, dt/t, j*xi*dt/(1 - j*xi*t)

# Let's check what the numeric value of each word is
# First, let's load the values from the pickle file

import pickle
vals_path = 'research/08_theory/A3f6_level16/vals8_150digit.pickle'
with open(vals_path, 'rb') as f:
    vals = pickle.load(f)

print("Type of vals:", type(vals))
print("Number of words in vals:", len(vals))
print("\nFirst few entries:")
for i, (k, v) in enumerate(list(vals.items())[:10]):
    print(f"  {k}: {v}")