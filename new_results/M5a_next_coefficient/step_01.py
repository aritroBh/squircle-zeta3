import os
# Explore the work directory structure
base = '<project>'
for root, dirs, files in os.walk(base):
    depth = root.replace(base, '').count(os.sep)
    if depth <= 2:
        print(root)
        if depth <= 1:
            for d in sorted(dirs):
                if d not in ['.git', '.venv', '__pycache__', 'node_modules']:
                    pass
