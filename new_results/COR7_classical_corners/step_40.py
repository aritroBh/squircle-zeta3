import os
base = '<project>/research/08_theory/COR5_corner_classicality/'
hits = []
for f in sorted(os.listdir(base)):
    if f.endswith('.py'):
        src = open(base+f).read()
        if 'corner' in src or 'quad' in src:
            hits.append(f)
print(hits)
for f in hits[:6]:
    print('='*20, f, '='*20)
    print(open(base+f).read()[:2500])
