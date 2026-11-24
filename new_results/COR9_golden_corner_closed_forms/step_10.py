base = '<project>/research/'
G1 = base+'10_frontier/G1_zeta3_any_angle/'
import json
for f in ['hpl_fit_k3_(0,0,0).json','hpl_fit_k3_(1,0,0).json','hpl_fit_k3_(0,1,0).json','hpl_fit_k2_(0,0,0).json','hpl_fit_k2_(1,0,0).json','hpl_fit_k1_(0,0,0).json','hpl_fit_k0_(0,0,0).json']:
    d = json.load(open(G1+f))
    print(f, len(d), d[:6])
    print()
