import os
base = '<project>'
for root, dirs, files in os.walk(os.path.join(base, 'research')):
    depth = root.replace(base, '').count(os.sep)
    if depth <= 3:
        print('  '*depth + os.path.basename(root) + '/')
        if depth <= 2:
            for f in files:
                print('  '*(depth+1) + f)
