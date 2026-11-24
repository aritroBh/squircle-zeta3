base = '<project>/research/'
print(open(base+'08_theory/COR4_M3alpha_classical/hpl_series.py').read())
import json
print(json.dumps(json.load(open(base+'10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(1,0,0).json')), indent=0))
print(json.dumps(json.load(open(base+'10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,1,0).json')), indent=0))
