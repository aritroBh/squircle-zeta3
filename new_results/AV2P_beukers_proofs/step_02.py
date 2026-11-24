import os
# Find AV2 directory
for root, dirs, files in os.walk('..'):
    if 'AV2' in root:
        print(root)
        for f in files:
            print('   ', f)
