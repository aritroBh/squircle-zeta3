import os
base = '<project>'
for d in ['research/08_theory/ZLAD1_zeta_ladder',
          'research/10_frontier/G1_zeta3_any_angle',
          'research/08_theory/G2_M4alpha',
          'research/08_theory/LEAD1_leading_genfun']:
    p = os.path.join(base, d)
    print('='*80)
    print(p, os.path.exists(p))
    if os.path.exists(p):
        for root, dirs, files in os.walk(p):
            print(root)
            for f in sorted(files):
                full = os.path.join(root, f)
                print('   ', f, os.path.getsize(full))
