import os
base = '<project>/research/08_theory'
k1e = os.path.join(base, 'K1e_hexagon_closed')
with open(os.path.join(k1e, 'step_11.py')) as f:
    c = f.read()
# print lines that mention lA or log or sub
for i, ln in enumerate(c.split('\n')):
    if any(s in ln for s in ['lA0','lA1','lB','log','subs','geo','w0','w1','wo','LA']):
        print(f"{i:3d}: {ln.rstrip()[:160]}")
