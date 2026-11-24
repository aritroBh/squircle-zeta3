import os
base = '<project>/research/08_theory'
e5c = os.path.join(base, 'E5c_nonlocal_z5')
with open(os.path.join(e5c, 'result.md')) as f:
    print(f.read())
