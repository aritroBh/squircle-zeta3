import os
tow2 = '<project>/research/08_theory/TOW2_tower_all_depths'
for f in sorted(os.listdir(tow2)):
    print(f, os.path.getsize(os.path.join(tow2,f)))
