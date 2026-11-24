import os
base = '<project>/research/08_theory'
with open(os.path.join(base,'HEXC2R_referee_certificate','result.md')) as f:
    print(f.read())
