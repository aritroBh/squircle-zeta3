import os
base = '<project>/research/08_theory'
e5e = os.path.join(base, 'E5e_pslq_340')
with open(os.path.join(e5e, 'result.md')) as f:
    print(f.read())
