import os
k1f = '<project>/research/08_theory/K1f_hexagon_order5'
with open(os.path.join(k1f, 'result.md')) as f:
    content = f.read()
print(content)
