G2='<project>/research/08_theory/G2_M4alpha/'
import json
d = json.load(open(G2+'hpl_fit_k4_(0,0,0).json'))
print("type", type(d), "len", len(d))
print("first 6 entries:", d[:6])
d2 = json.load(open(G2+'hpl_fit_k4_(1,0,0).json'))
print("psi2 first 6:", d2[:6])