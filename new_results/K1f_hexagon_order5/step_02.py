
import os
base = '<project>/research/08_theory'
for d in ['K1f_hexagon_order5', 'K1e_hexagon_closed', 'K1d_engine_fixed', 'K1c_engine_predict', 'E5_edge_law', 'E5c_nonlocal_z5']:
    p = os.path.join(base, d)
    print(f'=== {d} ===')
    for root, dirs, files in os.walk(p):
        dirs[:] = [x for x in dirs if x != '__pycache__']
        depth = root.replace(p, '').count(os.sep)
        indent = ' ' * 2 * depth
        print(f'{indent}{os.path.basename(root)}/')
        subindent = ' ' * 2 * (depth + 1)
        for f in sorted(files):
            print(f'{subindent}{f}')
    print()
