import os
p = '<project>/research/08_theory/A3h4_frame_ds'
for f in sorted(os.listdir(p)):
    if not f.startswith('.'): print(f, os.path.getsize(os.path.join(p,f)))
