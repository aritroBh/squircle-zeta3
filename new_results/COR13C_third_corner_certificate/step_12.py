import json
D = json.load(open('../../10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'))
print("k3 (0,0,0) fit: number of terms:", len(D))
for k, v in D: print("  ", k, "=", v)
