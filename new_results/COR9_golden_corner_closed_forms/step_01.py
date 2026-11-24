import os, glob
base = '<project>/research/'
for d in ['08_theory/COR7_classical_corners', '08_theory/COR5_corner_classicality', '10_frontier/G1_zeta3_any_angle', '08_theory/COR9_golden_corner_closed_forms']:
    p = os.path.join(base, d)
    print('='*80)
    print(p, os.path.isdir(p))
    if os.path.isdir(p):
        for f in sorted(glob.glob(p+'/*')):
            print('   ', os.path.basename(f), os.path.getsize(f))
