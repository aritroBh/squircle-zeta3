import os
# Explore cwd
for root, dirs, files in os.walk('.'):
    depth = root.count(os.sep)
    if depth > 2: 
        dirs[:] = []
        continue
    print(root, '->', files[:20])
