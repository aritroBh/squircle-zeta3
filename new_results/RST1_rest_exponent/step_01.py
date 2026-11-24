import os
os.chdir('<project>')
for root, dirs, files in os.walk('research/08_theory'):
    # limit depth
    depth = root.count(os.sep)
    if depth <= 3:
        print(root)
        for f in files:
            print('   ', f)
