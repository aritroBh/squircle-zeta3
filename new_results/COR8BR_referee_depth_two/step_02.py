import os
target = '<project>/research/08_theory/COR8B_depth_two_two_terms'
for root, dirs, files in os.walk(target):
    level = root.replace(target, '').count(os.sep)
    indent = ' ' * 2 * level
    print(f'{indent}{os.path.basename(root)}/')
    subindent = ' ' * 2 * (level + 1)
    for f in files:
        print(f'{subindent}{f}')
