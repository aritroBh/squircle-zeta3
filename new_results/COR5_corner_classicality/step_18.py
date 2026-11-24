import os, json
base = "<project>/research"
g1 = os.path.join(base, "10_frontier/G1_zeta3_any_angle")
with open(os.path.join(g1, "values_2pi_3.txt")) as f:
    print(f.read()[:3000])
