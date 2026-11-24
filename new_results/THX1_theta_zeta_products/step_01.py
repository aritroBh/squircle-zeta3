
import os, subprocess
print(os.getcwd())
for root, dirs, files in os.walk('.'):
    # skip deep
    depth = root.count(os.sep)
    if depth > 2:
        dirs[:] = []
        continue
    for f in files:
        p = os.path.join(root, f)
        print(p, os.path.getsize(p))
