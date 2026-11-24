import os
base = '<project>/research/08_theory'
for root, dirs, files in os.walk(base):
    # only show relevant dirs
    if 'WT2' in root:
        print(root)
        for f in files:
            print('   ', f)
