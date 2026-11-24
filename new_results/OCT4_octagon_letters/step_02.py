import os
base = 'research/08_theory/'
for d in ['OCT4_octagon_letters', 'D15b_polygon_law', 'OCT3_octagon_basis']:
    p = os.path.join(base, d)
    print('='*60)
    print(d, '->', os.path.exists(p))
    if os.path.exists(p):
        for f in sorted(os.listdir(p)):
            full = os.path.join(p, f)
            sz = os.path.getsize(full) if os.path.isfile(full) else '<dir>'
            print('   ', f, sz)
