
import os
base='research/08_theory'
k1d=os.path.join(base,'K1d_engine_fixed')
with open(os.path.join(k1d,'esum_engine.py')) as f:
    src=f.read()
print(len(src))
print(src[:5500])
