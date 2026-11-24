import os, subprocess
# Explore the working directory and files under review
base = os.getcwd()
print("CWD:", base)
for root, dirs, files in os.walk(base):
    # avoid deep recursion into unrelated dirs
    level = root.replace(base, '').count(os.sep)
    indent = ' ' * 2 * level
    print(f'{indent}{os.path.basename(root)}/')
    subindent = ' ' * 2 * (level + 1)
    for f in files:
        print(f'{subindent}{f}')
    if level > 3:
        dirs[:] = []
