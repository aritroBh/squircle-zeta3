import os
base = "<project>/research/08_theory"
for root, dirs, files in os.walk(base):
    for f in files:
        print(os.path.join(root, f))
