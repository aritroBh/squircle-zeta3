import os
for root, dirs, files in os.walk('.'):
    # skip .git
    if '.git' in root: continue
    depth = root.count(os.sep)
    if depth > 2: continue
    for f in files:
        print(os.path.join(root, f))
