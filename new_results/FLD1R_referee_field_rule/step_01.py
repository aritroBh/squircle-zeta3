import os
print(os.getcwd())
for root, dirs, files in os.walk('.'):
    depth = root.count(os.sep)
    if depth < 3:
        print(root, '->', files[:20])
