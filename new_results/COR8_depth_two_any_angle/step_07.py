
import json, re
base='<project>/research/'
G1=base+'10_frontier/G1_zeta3_any_angle/'
R3=json.load(open(G1+'hpl_fit_k3_(0,0,0).json'))
print(len(R3))
for n,c in R3: print(n, c)
