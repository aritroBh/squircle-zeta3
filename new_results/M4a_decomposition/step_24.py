import os
BASE = '<project>'
OUT = os.path.join(BASE,'paper/squircle_note/outputs')
for f in ['mk_numerics_A_1200.out','mk_numerics_C_60.out','m4_tails_exact.pkl']:
    fp = os.path.join(OUT,f)
    if f.endswith('.out'):
        print('='*20, f)
        print(open(fp).read())
import pickle
d = pickle.load(open(os.path.join(OUT,'m4_tails_exact.pkl'),'rb'))
print('m4_tails_exact.pkl keys:', list(d.keys()))
for k,v in d.items():
    print(k, type(v))
