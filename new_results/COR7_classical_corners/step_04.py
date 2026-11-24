import os
base = '<project>/research/'
for d in ['08_theory/COR6_true_point_symbols/', '08_theory/COR5_corner_classicality/', '10_frontier/G1_zeta3_any_angle/']:
    root = base + d
    print(root)
    for f in sorted(os.listdir(root)):
        p = os.path.join(root, f)
        print('   ', f, os.path.getsize(p))
    print()
