import os, subprocess
base = '<project>'
# find ZQ5 dirs
for root, dirs, files in os.walk(base):
    if '.venv' in root: continue
    for d in list(dirs):
        if 'ZQ5' in d:
            print(os.path.join(root,d))
    # limit depth
    if root.count(os.sep) - base.count(os.sep) > 3:
        dirs[:] = []
