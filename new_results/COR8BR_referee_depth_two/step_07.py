import os
ref = '<project>/research/08_theory/COR8R_referee_depth_two'
with open(os.path.join(ref, 'result.md')) as f:
    print(f.read())
