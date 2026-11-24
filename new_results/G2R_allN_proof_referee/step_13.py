import json
base = "<project>/research/08_theory/G2_M4alpha/"
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
for f in ["hpl_fit_k4_(0,0,0).json", "hpl_fit_k4_(1,0,0).json"]:
    d = json.load(open(base+f))
    print(f, len(d), "terms; first 5:", d[:5])
import os
print(os.path.exists(G1+"hpl_fit_k3_(1,0,0).json"), os.path.exists(G1+"flat_exact_N320_K3.json"))
d = json.load(open(G1+"hpl_fit_k3_(1,0,0).json"))
print("G1 k3 (1,0,0):", len(d), d[:5])
# check flat_exact_N640_K4_formula.json structure
d = json.load(open(base+"flat_exact_N640_K4_formula.json"))
print(d.keys(), list(d['4'].items())[:2])
