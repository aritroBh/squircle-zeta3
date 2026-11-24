import os
base='<project>/research/08_theory'
src=open(os.path.join(base,'WT2_weight2_all_k/fastred.py')).read()
print(len(src))
# find functions pf_sub, pf_scale, pf_check_reality
import re
for m in re.finditer(r'def (\w+)', src):
    print(m.group(1), m.start())
