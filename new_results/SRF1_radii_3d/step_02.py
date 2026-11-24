
import os
d16 = '<project>/research/08_theory'
# find D16 related files
for root, dirs, files in os.walk(d16):
    for f in files:
        if 'D16' in f or 'd16' in f or 'pure_weight' in root:
            print(os.path.join(root, f))
