import os, subprocess
print(os.getcwd())
print(os.listdir('.'))
# check the paper location
for root, dirs, files in os.walk('../../..', topdown=True):
    depth = root.count(os.sep)
    if depth < 3:
        print(root, dirs[:8], files[:8])
    if depth > 2:
        dirs[:] = []
