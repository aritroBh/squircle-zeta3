
base = "<project>/research/08_theory/"
import os
p = base+"DEC5D_decagon_order5_depth.brief"
print(os.path.exists(p))
if os.path.exists(p):
    print(open(p).read()[:8000])
