import os
base = '<project>/research/08_theory'
k1e = os.path.join(base, 'K1e_hexagon_closed')
with open(os.path.join(k1e, 'result.md')) as f:
    print(f.read())
