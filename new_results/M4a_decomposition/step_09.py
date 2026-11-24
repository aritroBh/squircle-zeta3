import os
BASE = '<project>'
OUT = os.path.join(BASE,'paper/squircle_note/outputs')
for f in ['ibp_reduce_k4.out','m4_tails_exact.out','ibp_reduce_k123.out']:
    print('='*30, f)
    print(open(os.path.join(OUT,f)).read())
