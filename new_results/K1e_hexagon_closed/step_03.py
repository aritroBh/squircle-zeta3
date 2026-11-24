
import os
base='research/08_theory'
k1e=os.path.join(base,'K1e_hexagon_closed')
print(os.listdir(k1e))
for f in ['K1e_hexagon_closed.brief']:
    with open(os.path.join(base,f)) as fh:
        print(fh.read()[:3000])
