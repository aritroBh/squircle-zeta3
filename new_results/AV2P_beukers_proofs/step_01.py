import os
# Explore directory structure
for root, dirs, files in os.walk('..'):
    level = root.replace('..', '').count(os.sep)
    if level <= 2:
        indent = ' ' * 2 * level
        print(f'{indent}{os.path.basename(root)}/')
        subindent = ' ' * 2 * (level + 1)
        for f in files:
            print(f'{subindent}{f}')
