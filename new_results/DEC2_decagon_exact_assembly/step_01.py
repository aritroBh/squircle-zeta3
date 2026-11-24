import os
base = "<project>"
for root, dirs, files in os.walk(os.path.join(base, "research/08_theory")):
    depth = root.replace(base, "").count(os.sep)
    if depth <= 4:
        print(root)
        for f in files:
            print("   ", f)
