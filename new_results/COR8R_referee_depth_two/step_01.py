import os, subprocess
# Find the base research directory
cwd = os.getcwd()
print("CWD:", cwd)
base = cwd
# go up until we find 'research'
while not os.path.isdir(os.path.join(base, 'research')) and base != '/':
    base = os.path.dirname(base)
print("BASE:", base)
for root, dirs, files in os.walk(os.path.join(base, 'research', '08_theory')):
    depth = root.replace(os.path.join(base,'research','08_theory'), '').count(os.sep)
    if depth <= 2:
        print(root)
        for f in files:
            print('   ', f)
