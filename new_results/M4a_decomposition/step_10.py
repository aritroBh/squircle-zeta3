import os
BASE = '<project>'
S = os.path.join(BASE,'paper/squircle_note/scripts')
for f in ['ibp_reduce.py','mk_integrands.py','qsqrt2.py']:
    print('='*30, f)
    print(open(os.path.join(S,f)).read())
