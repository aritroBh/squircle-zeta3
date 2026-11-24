import os
base = '<project>/research/08_theory'
e5 = os.path.join(base, 'E5_edge_law')
with open(os.path.join(e5, 'result.md')) as f:
    print(f.read())
