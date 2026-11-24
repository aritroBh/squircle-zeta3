import os, shutil, pickle, sys

base = '<home>)'
a3f3 = base + '/research/08_theory/A3f3_level8_compute'
print(os.path.isdir(a3f3))
for root, dirs, files in os.walk(a3f3):
    depth = root[len(a3f3):].count(os.sep)
    if depth > 1:
        continue
    print(root, files[:20])
