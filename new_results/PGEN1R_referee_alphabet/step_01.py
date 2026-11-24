import os, subprocess
print(os.getcwd())
for root, dirs, files in os.walk('.'):
    if '.git' in root: continue
    depth = root.count(os.sep)
    if depth < 2:
        print(root, files[:20])
