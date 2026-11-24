base = '<project>/research/08_theory/'
import os
# CX8 review
for cand in ['CX8_tow2_review/report.md']:
    p = base+cand
    print(p, os.path.exists(p))
# find CX8 directory
cx = [d for d in os.listdir(base) if d.startswith('CX8')]
print(cx)
