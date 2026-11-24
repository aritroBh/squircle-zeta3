
import os
base = "<project>"
p = os.path.join(base,"research/08_theory/SRF3R_referee_octahedron/result.md")
print(os.path.exists(p))
print(open(p).read())
