import os
base = os.getcwd()
print(base)
for root, dirs, files in os.walk('.'):
    depth = root.count(os.sep)
    if depth <= 2:
        print(root, '->', files[:20])
