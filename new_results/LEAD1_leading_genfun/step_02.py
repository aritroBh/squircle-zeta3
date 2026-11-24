import os
base = '<project>'
wts = os.path.join(base, 'research', '08_theory', 'WTS_weight_k_structure')
for root, dirs, files in os.walk(wts):
    depth = root.replace(wts, '').count(os.sep)
    if depth <= 2:
        print(root)
        for f in sorted(files)[:40]:
            print('   ', f)
