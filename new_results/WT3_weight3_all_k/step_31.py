
import os, pickle, sys
BASE = '<project>'
WT2B = os.path.join(BASE, 'research/08_theory/WT2b_stabilization')
with open(os.path.join(WT2B,'record_W2b.pkl'),'rb') as fh:
    rec = pickle.load(fh)
print(type(rec))
if isinstance(rec, dict):
    for kk in rec:
        v = rec[kk]
        print(kk, type(v), (len(v) if hasattr(v,'__len__') else v))
