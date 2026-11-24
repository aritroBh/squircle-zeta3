import os, json
base = "<project>/research"
g1 = os.path.join(base, "10_frontier/G1_zeta3_any_angle")
R3 = json.load(open(os.path.join(g1, "hpl_fit_k3_(0,0,0).json")))
for n, c in R3:
    print(f"{c:>8}  {n}")
