import os
print(os.getcwd())
for root, dirs, files in os.walk('.'):
    # skip hidden dirs
    dirs[:] = [d for d in dirs if not d.startswith('.')]
    level = root.count(os.sep)
    if level > 3: 
        dirs[:] = []
        continue
    indent = ' ' * 2 * level
    print(f'{indent}{os.path.basename(root)}/')
    subindent = ' ' * 2 * (level + 1)
    for f in files:
        if not f.startswith('.'):
            print(f'{subindent}{f}')
