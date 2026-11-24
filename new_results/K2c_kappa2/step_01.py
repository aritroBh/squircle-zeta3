import os
base = '<project>/research/08_theory'
for root, dirs, files in os.walk(base):
    depth = root.replace(base,'').count(os.sep)
    if depth <= 2:
        print(root)
        for f in files:
            if f.endswith(('.md','.py')):
                print('   ', f)
