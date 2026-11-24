import json, os
base = '<project>'
p = os.path.join(base,'research/10_frontier/G1_zeta3_any_angle')
print(sorted(os.listdir(p)))
data = json.load(open(os.path.join(p,'hpl_fit_k3_(0,0,0).json')))
print(type(data))
if isinstance(data, dict):
    print(list(data.keys())[:20])
