import os
for p in ['../../10_frontier/F1_fourier_squircle', '../../10_frontier/H1_octahedron_end', '../D11c_zeta_first_order', '../D11_rhombus_resurgence']:
    ap = os.path.abspath(p)
    print('====', ap, os.path.exists(ap))
    if os.path.exists(ap):
        for root, dirs, files in os.walk(ap):
            print(root, files[:30])
