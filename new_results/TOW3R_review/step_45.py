import os, subprocess
base = '<project>/research/08_theory/'
# find TOW1 proof docs
for d in ['TOW1_tower_law']:
    print(os.listdir(base+d))
