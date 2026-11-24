import os, glob

base = "<project>"
# Read SRF1 result.md
with open(os.path.join(base, "research/08_theory/SRF1_radii_3d/result.md")) as f:
    print(f.read())
