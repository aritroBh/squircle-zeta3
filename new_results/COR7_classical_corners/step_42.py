import os
base = '<project>/research/08_theory/COR5_corner_classicality/'
hits = []
for f in sorted(os.listdir(base)):
    if f.endswith('.py'):
        src = open(base+f).read()
        if '1100' in src or 'expm' in src or 'extraprec' in src or 'small' in src.lower() or 'gauss' in src.lower():
            hits.append((f, [l for l in src.split('\n') if '1100' in l or 'extraprec' in l or 'gauss' in l.lower()][:3]))
for f, h in hits:
    print(f, h)
