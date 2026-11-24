import os
base = os.path.abspath('..')
for d in sorted(os.listdir(base)):
    print(d)
