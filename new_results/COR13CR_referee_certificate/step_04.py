import os
base = '<project>/research/08_theory'
with open(os.path.join(base,'COR13C_third_corner_certificate','result.md')) as f:
    print(f.read())
