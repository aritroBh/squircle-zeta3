import os, pickle
base = os.path.dirname(os.getcwd())
d2 = os.path.join(base,'DOD2_dodecagon_finish')
d1 = os.path.join(base,'DOD1_dodecagon')

def peek(path):
    with open(path,'rb') as f:
        obj = pickle.load(f)
    print('---', os.path.basename(path), type(obj))
    if isinstance(obj, dict):
        for k,v in list(obj.items())[:12]:
            print('   key:', repr(k)[:80], '->', type(v), repr(v)[:100])
    elif isinstance(obj, (list,tuple)):
        for v in list(obj)[:8]:
            print('   ', type(v), repr(v)[:100])
    else:
        print(repr(obj)[:300])

for f in ['Xi_verified_620.pkl','atom_coefs_pq.pkl','gdswaps620.pkl']:
    peek(os.path.join(d2,f))
for f in ['P3_630.pkl','S_pq_620.pkl','gamdel620.pkl','atoms620.pkl']:
    peek(os.path.join(d1,f))
