import os, json

# Explore the data files
base = 'research/08_theory/'
for root, dirs, files in os.walk(base):
    # limit depth
    depth = root.replace(base, '').count(os.sep)
    if depth <= 2:
        print(root)
        for f in files:
            print('   ', f)
