
import os, subprocess
base = "research/08_theory"
for root, dirs, files in os.walk(base):
    # limit depth
    depth = root.count(os.sep)
    if depth <= 3:
        print(root, "->", files[:20])
