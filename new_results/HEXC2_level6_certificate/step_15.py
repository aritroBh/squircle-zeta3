import json
from collections import Counter
G1 = '<project>/research/10_frontier/G1_zeta3_any_angle'
def load_fit(name):
    D = json.load(open(G1+'/'+name))
    out=[]
    for key,val in D:
        m = re.match(r"c\^(-?\d+) H\[(.*)\]", key)
        if key=="1": e,w = 0,()
        elif key=="c": e,w = 1,()
        else:
            e=int(m.group(1)); w=tuple(int(z) for z in m.group(2).split(',')) if m.group(2).strip() else ()
        out.append((e,w,val))
    return out
import re
R3 = load_fit('hpl_fit_k3_(0,0,0).json')
PSI = load_fit('hpl_fit_k3_(1,0,0).json')
by = {}
for e,w,v in R3: by.setdefault(len(w),[]).append((e,w,v))
for k in sorted(by):
    print('R3 weight',k,'terms:',len(by[k]))
    for e,w,v in by[k]: print('   c^%d H%s  %s'%(e,list(w),v))
