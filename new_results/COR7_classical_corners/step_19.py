import json
base = '<project>/research/'
D = json.load(open(base+'10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'))
print(len(D), "terms in R3")
for k, v in D:
    print(k, v)
