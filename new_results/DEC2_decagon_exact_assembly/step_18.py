base = "<project>/research/08_theory"
# find the actual builder file step_40 in K1d
import os
p = f"{base}/K1d_engine_fixed"
# step_40.py printed above was actually the result-writer. Look for the builder
with open(f"{p}/step_33_builder.py") as f:
    print(f.read()[:3000])
