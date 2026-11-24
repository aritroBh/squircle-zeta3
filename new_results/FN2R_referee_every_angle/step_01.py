import os, subprocess
base = "research/08_theory"
for root, dirs, files in os.walk(base):
    depth = root.replace(base, "").count(os.sep)
    if depth <= 2:
        print(root, "->", files[:20])
