
import os
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
print(os.path.exists(G1 + "hpl_fit_k3_(1,0,0).json"))
import json
d = json.load(open(G1 + "hpl_fit_k3_(1,0,0).json"))
print(len(d), "terms;", d[:6])
