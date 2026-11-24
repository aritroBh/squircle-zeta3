
import os
chn1 = '<project>/research/08_theory/CHN1_chain_bodies'
for f in sorted(os.listdir(chn1)):
    p = os.path.join(chn1, f)
    sz = os.path.getsize(p)
    print(f"  {f}   {sz}")
