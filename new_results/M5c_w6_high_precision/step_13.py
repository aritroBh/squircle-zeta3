import os
base='<project>/research/08_theory/M5b_w6_closed_form'
log=open(os.path.join(base,'progress.log')).read().splitlines()
import re
for l in log:
    if re.search(r'w6 (bits|complete)|lchi|Lchi|L\(χ|seconds=|integral', l, re.I):
        print(l)
