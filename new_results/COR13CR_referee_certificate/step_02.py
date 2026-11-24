import os
base = '<project>/research/08_theory'
for d in ['COR13C_third_corner_certificate', 'HEXC2R_referee_certificate']:
    p = os.path.join(base, d)
    print("="*80)
    print(d, os.path.exists(p))
    if os.path.exists(p):
        for f in sorted(os.listdir(p)):
            print("  ", f, os.path.getsize(os.path.join(p,f)))
