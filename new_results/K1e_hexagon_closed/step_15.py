
import os
base='research/08_theory'
k1c=os.path.join(base,'K1c_engine_predict')
with open(os.path.join(k1c,'step_01.py')) as f:
    print(f.read()[:5000])
