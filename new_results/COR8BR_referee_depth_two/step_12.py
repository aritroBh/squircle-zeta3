import json, subprocess
p = '<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'
with open(p) as f:
    data = json.load(f)
print(type(data))
if isinstance(data, dict):
    print(list(data.keys()))
    for k in list(data.keys())[:5]:
        print(k, str(data[k])[:500])
elif isinstance(data, list):
    print(len(data))
    for x in data[:5]:
        print(str(x)[:300])
