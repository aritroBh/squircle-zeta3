import json, os
base = '<project>'
data = json.load(open(os.path.join(base,'research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json')))
print("num entries:", len(data))
print(json.dumps(data[:6], indent=1)[:3000])
