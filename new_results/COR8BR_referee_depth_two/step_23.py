import os
c8 = '<project>/research/08_theory/COR8_depth_two_any_angle'
print(os.listdir(c8))
for f in os.listdir(c8):
    if f.endswith('.md'):
        with open(os.path.join(c8, f)) as fh:
            print(fh.read()[:8000])
