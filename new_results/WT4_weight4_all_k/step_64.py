import os
work='<project>/research/08_theory/WT4_weight4_all_k'
src=open(os.path.join(work,'verify_wt4.py')).read()
# find remaining top-level uses of aN/bN/r2 outside quad_value
import re
for i, line in enumerate(src.splitlines()):
    if re.search(r'\b(aN|bN)\b', line) or re.search(r'\br2\b', line):
        print(i, line[:110])
