base = "<project>/research/08_theory"
import os
p = f"{base}/K1d_engine_fixed"
print(sorted(os.listdir(p)))
with open(f"{p}/step_40.py") as f:
    print(f.read())
