
import os, pickle, sys
BASE = '<project>'
WT2B = os.path.join(BASE, 'research/08_theory/WT2b_stabilization')
rec = pickle.load(open(os.path.join(WT2B,'record_W2b.pkl'),'rb'))
e = rec[(4,(3,0))]
print(e.keys())
for kk, vv in e.items():
    print("---", kk, type(vv))
    print(vv)
