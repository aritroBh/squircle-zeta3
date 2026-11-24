import os
base = '<project>'
for root, dirs, files in os.walk(base):
    if 'PRV' in root or any('PRV' in f for f in files):
        print(root)
        for f in sorted(files):
            print('   ', f)
