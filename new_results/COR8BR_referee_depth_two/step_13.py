import json
p = '<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'
with open(p) as f:
    data = json.load(f)
for x in data:
    print(x)
