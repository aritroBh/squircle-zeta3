import os
# Explore working directory structure
for root, dirs, files in os.walk('..'):
    depth = root.count(os.sep)
    if depth <= 2:
        print(root)
        for f in files:
            if f.endswith(('.py','.md','.json','.txt')):
                print('   ', f)
