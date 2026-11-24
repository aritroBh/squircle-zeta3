
import os, subprocess
# First, find my work directory and explore the research tree
base = '<project>'
for root, dirs, files in os.walk(base):
    # skip hidden and venv
    dirs[:] = [d for d in dirs if not d.startswith('.') and d != '__pycache__']
    depth = root.replace(base, '').count(os.sep)
    if depth <= 3:
        indent = ' ' * 2 * depth
        print(f'{indent}{os.path.basename(root)}/')
        if depth <= 2:
            subindent = ' ' * 2 * (depth + 1)
            for f in files[:20]:
                print(f'{subindent}{f}')
