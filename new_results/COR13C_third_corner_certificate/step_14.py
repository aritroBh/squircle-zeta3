import json
for f in ['hpl_fit_k3_(1,0,0).json', 'hpl_fit_k2_(1,0,0).json', 'hpl_fit_k3_(0,1,0).json']:
    D = json.load(open('../../10_frontier/G1_zeta3_any_angle/' + f))
    print(f, ":", len(D), "terms")
    for k, v in D: print("   ", k, "=", v)
    print()
