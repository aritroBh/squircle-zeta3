import os, glob
base = '<project>'
# explore structure
for root, dirs, files in os.walk(os.path.join(base, 'research', '08_theory')):
    depth = root.replace(os.path.join(base,'research','08_theory'), '').count(os.sep)
    if depth <= 2:
        print(root)
        for f in files[:20]:
            print('   ', f)
