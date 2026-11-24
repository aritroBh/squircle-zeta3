
import os
base = ".."
for root, dirs, files in os.walk(base):
    depth = root.count(os.sep)
    if depth <= 2:
        print(root, "->", sorted(files))
