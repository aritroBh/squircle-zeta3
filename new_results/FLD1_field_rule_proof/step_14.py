import json
base='<project>/research/10_frontier/G1_zeta3_any_angle/'
for fn in ['hpl_fit_k2_(0,0,0).json','hpl_fit_k3_(0,0,0).json','hpl_fit_k3_(1,0,0).json','hpl_fit_k3_(0,1,0).json']:
    print("="*70); print(fn)
    data = json.load(open(base+fn))
    print(f"  ({len(data)} terms)")
    for nm, v in data:
        print(f"   {v}  *  {nm}")