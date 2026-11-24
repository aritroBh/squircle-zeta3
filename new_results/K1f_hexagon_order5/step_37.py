import os, time
wd = '<project>/research/08_theory/K1f_hexagon_order5'
for f in ['step_24.out','step_24.log','E5_polys.pkl']:
    p = f'{wd}/{f}'
    if os.path.exists(p):
        print(f, os.path.getsize(p))
        if f.endswith('.out') or f.endswith('.log'):
            print(open(p).read())
    else:
        print(f, 'MISSING')
