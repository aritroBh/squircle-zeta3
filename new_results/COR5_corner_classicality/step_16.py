import os
base = "<project>/research"
g1 = os.path.join(base, "10_frontier/G1_zeta3_any_angle")
with open(os.path.join(g1, "m3_hpl_check.py")) as f:
    print(f.read())
