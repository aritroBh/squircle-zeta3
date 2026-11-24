
import os
base='research/08_theory'
k1e=os.path.join(base,'K1e_hexagon_closed')
for f in sorted(os.listdir(k1e)):
    print('='*30, f)
    with open(os.path.join(k1e,f)) as fh:
        print(fh.read()[:2500])
