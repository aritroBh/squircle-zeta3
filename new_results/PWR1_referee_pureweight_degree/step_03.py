import os
# Find the referenced result.md files
paths = [
    '../D16_pure_weight/result.md',
    '../../10_frontier/H1_octahedron_end/result.md',
    '../K3_fourier_all_orders/result.md',
]
for p in paths:
    ap = os.path.abspath(p)
    print(ap, os.path.exists(ap))
# find the paper
for root, dirs, files in os.walk('../../..'):
    if 'main.tex' in files and 'ieee' in root:
        print('PAPER:', os.path.join(root, 'main.tex'))
        break
