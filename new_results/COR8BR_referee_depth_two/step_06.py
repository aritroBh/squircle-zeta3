import os
ref = '<project>/research/08_theory/COR8R_referee_depth_two'
for root, dirs, files in os.walk(ref):
    level = root.replace(ref, '').count(os.sep)
    indent = ' ' * 2 * level
    print(f'{indent}{os.path.basename(root)}/')
    subindent = ' ' * 2 * (level + 1)
    for f in files:
        print(f'{subindent}{f}')
