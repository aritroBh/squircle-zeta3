import os
h4 = '<project>/research/08_theory/A3h4_frame_ds'
for f in sorted(os.listdir(h4)):
    print(f, os.path.getsize(os.path.join(h4,f)))
