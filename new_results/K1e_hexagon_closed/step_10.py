
import os
base='research/08_theory'
k1d=os.path.join(base,'K1d_engine_fixed')
with open(os.path.join(k1d,'step_42.py')) as f:
    print(f.read())
