import os, json
base = '<project>/research/08_theory'
d = os.path.join(base,'COR13C_third_corner_certificate')
for f in ['fast_cmzv.py','numeric_cmzv.py','exact_linalg.py']:
    print("="*30, f)
    print(open(os.path.join(d,f)).read()[:3500])
