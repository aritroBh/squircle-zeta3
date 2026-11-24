
import os
base='research/08_theory'
k1d=os.path.join(base,'K1d_engine_fixed')
with open(os.path.join(k1d,'result.md')) as f:
    txt=f.read()
print(len(txt))
print(txt[:6000])
