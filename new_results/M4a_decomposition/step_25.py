import os, sys, pickle
BASE = '<project>'
SCR = os.path.join(BASE,'paper/squircle_note/scripts')
sys.path.insert(0, SCR)
OUT = os.path.join(BASE,'paper/squircle_note/outputs')
d = pickle.load(open(os.path.join(OUT,'m4_tails_exact.pkl'),'rb'))
print('keys:', list(d.keys()))
for k,v in d.items():
    print(k, type(v), repr(v)[:200])
