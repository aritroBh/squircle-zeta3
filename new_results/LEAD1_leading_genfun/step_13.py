import os, pickle, json
base = '<project>'
theory = os.path.join(base, 'research', '08_theory', 'WTS_weight_k_structure', 'theory')
with open(os.path.join(theory, 'kern_10.pkl'),'rb') as f:
    K10 = pickle.load(f)
rec = K10[(10,0)]
print(type(rec))
print(rec if isinstance(rec,(list,tuple)) else list(rec.keys()) if isinstance(rec,dict) else rec)
import pprint
if isinstance(rec, dict):
    for k,v in rec.items():
        print(k, type(v), str(v)[:300])
