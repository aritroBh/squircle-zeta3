import os
print(os.getcwd())
for root, dirs, files in os.walk('.'):
    # limit depth
    depth = root.count(os.sep)
    if depth < 3:
        print(root)
        for f in files:
            print('   ', f)
