import os, subprocess
print("cwd:", os.getcwd())
for root, dirs, files in os.walk('.'):
    # limit depth
    depth = root.count(os.sep)
    if depth <= 2:
        print(root, files[:20])
