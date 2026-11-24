import os
base = '<project>'
for root, dirs, files in os.walk(base):
    if 'K1f' in root:
        print(root, '::', len(files), 'files')
