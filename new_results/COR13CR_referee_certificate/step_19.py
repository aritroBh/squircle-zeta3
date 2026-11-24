import os, json
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
print(sorted(os.listdir(G1))[:40])
f = json.load(open(G1 + "/hpl_fit_k3_(0,0,0).json"))
print("fit k3 (0,0,0): entries", len(f))
for e in f[:10]: print(e)
