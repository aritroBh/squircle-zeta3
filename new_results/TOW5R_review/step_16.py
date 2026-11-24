t1 = '<project>/research/08_theory/TOW1_tower_law'
import os
print(sorted(os.listdir(t1)))
print()
# search for S(t) definition in the proof docs
for f in sorted(os.listdir(t1)):
    if f.endswith('.md'):
        s = open(os.path.join(t1,f)).read()
        print('=====', f, len(s))
