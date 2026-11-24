import os
base = "<project>/research"
g1 = os.path.join(base, "10_frontier/G1_zeta3_any_angle")
with open(os.path.join(g1, "m2_2pi3_pslq.py")) as f:
    print(f.read())
