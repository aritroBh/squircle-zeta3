
import os, pickle, json
from fractions import Fraction as F
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
os.chdir(WD)
K5 = pickle.load(open('kern_5.pkl', 'rb'))
PAIR_TO_POLE = {(0,0):'0', (1,1):'a', (-1,-1):'-a', (-1,1):'b', (1,-1):'-b'}

def kern_pkl_to_omega(K, keys):
    out = {}
    for key in keys:
        rec = []
        for (p0,p1), (c0,c1,c2,c3) in K[key]['real'].items():
            assert c2==0 and c3==0, (key, p0, p1, c2, c3)
            if c0==0 and c1==0: continue
            rec.append([PAIR_TO_POLE[(int(p0),int(p1))], [str(c0), str(c1)]])
        rec.sort(key=lambda t: ['0','a','-a','b','-b'].index(t[0]))
        out[f"{key[0]},{key[1]}"] = rec
    return out

om_k5 = kern_pkl_to_omega(K5, [(5,0),(4,1),(3,2)])
om5 = json.load(open('omega5.json'))
for key in ["5,0","4,1","3,2"]:
    mine = sorted((t[0],t[1][0],t[1][1]) for t in om_k5[key])
    theirs = sorted((t[0],t[1][0],t[1][1]) for t in om5[key])
    print(key, "kern_5.pkl == omega5.json:", mine==theirs)
print(json.dumps(om_k5, indent=1))
