import os
t5 = '<project>/research/08_theory/TOW5_tower_depth10'
for f in sorted(os.listdir(t5)):
    p = os.path.join(t5, f)
    print(f, os.path.getsize(p) if os.path.isfile(p) else '<dir>')
