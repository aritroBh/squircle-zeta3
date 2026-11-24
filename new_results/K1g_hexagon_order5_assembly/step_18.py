import os
base = '<project>/research/08_theory'
k1c = os.path.join(base, 'K1c_engine_predict')
with open(os.path.join(k1c, 'result.md')) as f:
    print(f.read())
