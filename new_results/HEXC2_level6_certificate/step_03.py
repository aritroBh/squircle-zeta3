import os
for d in ['HEXC1_hexagon_certificate', 'A3h5_exact_certificate', 'A3h4_frame_ds', 'COR5_corner_classicality']:
    p = os.path.join('<project>/research/08_theory', d)
    print('===', d)
    for f in sorted(os.listdir(p)):
        if not f.startswith('.'):
            full = os.path.join(p, f)
            print(f'  {f}  {os.path.getsize(full)} bytes')
