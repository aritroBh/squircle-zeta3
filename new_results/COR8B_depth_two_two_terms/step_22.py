
import json
d = json.load(open('../../10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'))
print(type(d))
if isinstance(d, dict):
    print(list(d.keys())[:20])
    import json as J
    s = J.dumps(d)[:3000]
    print(s)
elif isinstance(d, list):
    print("len", len(d))
    print(json.dumps(d[:5], indent=1)[:3000])
