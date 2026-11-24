import os
print(os.getcwd())
for root, dirs, files in os.walk('..'):
    # limit depth
    depth = root.count(os.sep)
    if depth <= 2:
        for f in files:
            print(os.path.join(root, f))
