import json, re
G1 = '<project>/research/10_frontier/G1_zeta3_any_angle'
import os
print(sorted(os.listdir(G1)))
R3 = json.load(open(G1+'/hpl_fit_k3_(0,0,0).json'))
PSI = json.load(open(G1+'/hpl_fit_k3_(1,0,0).json'))
print('R3 terms:', len(R3), ' PSI terms:', len(PSI))
for key, val in R3[:8]:
    print(key, val)
