import os
for name in ['OCT5_octagon_exact','K1d_engine_fixed','D15b_polygon_law','E1_edge_law','E5c_nonlocal_z5','E5_edge_law','D15c_edge_zeta3']:
    d=f'research/08_theory/{name}'
    print(name, ':', sorted(os.listdir(d)) if os.path.isdir(d) else 'MISSING')
