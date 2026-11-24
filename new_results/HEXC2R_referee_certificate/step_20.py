
import json
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
D = json.load(open(G1 + "/hpl_fit_k3_(0,0,0).json"))
print("R3 fit terms:", len(D))
for k, v in D[:12]:
    print(k, "=", v)
print("...")
D1 = json.load(open(G1 + "/hpl_fit_k3_(1,0,0).json"))
print("Psi fit terms:", len(D1))
for k, v in D1[:12]:
    print(k, "=", v)
