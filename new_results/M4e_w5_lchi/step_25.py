import os
p = '<project>/research/08_theory/M4b2_w4pp_nice'
print(sorted(os.listdir(p)))
if os.path.exists(os.path.join(p,'nice_display.md')):
    print(open(os.path.join(p,'nice_display.md')).read()[:5000])
