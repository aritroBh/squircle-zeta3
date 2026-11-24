import os
base = "<project>/research/10_frontier/G1_zeta3_any_angle/"
files = os.listdir(base)
print([f for f in files if 'm3' in f.lower() or 'hpl_fit_k3' in f])
# read result.md
with open(base+"result.md") as f:
    txt = f.read()
print(len(txt))
# find section B3
import re
idx = txt.find("B3")
print(txt[idx-200:idx+4000])
