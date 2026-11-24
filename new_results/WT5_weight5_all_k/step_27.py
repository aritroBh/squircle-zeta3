base='<project>'
import re
with open(base+'/research/08_theory/M5b_w6_closed_form/progress.log') as f:
    lines=f.read().splitlines()
for l in lines:
    if 'w5pp' in l or 'W5' in l or 'LLL' in l or 'symbol solve' in l or 'HNF' in l:
        print(l[:220])
