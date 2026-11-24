
import os
base = "<project>/research/08_theory/"
for d in sorted(os.listdir(base)):
    if 'EC' in d or 'dec' in d.lower():
        print(d, os.path.isdir(base+d) and sorted(os.listdir(base+d))[:40])
