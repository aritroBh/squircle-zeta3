import os
print(os.getcwd())
for f in sorted(os.listdir('.')):
    p = os.path.join('.',f)
    print(f, os.path.getsize(p))
