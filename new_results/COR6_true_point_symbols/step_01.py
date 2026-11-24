import os
base = "<project>/research/"
for root, dirs, files in os.walk(base + "08_theory"):
    depth = root.replace(base + "08_theory", "").count(os.sep)
    if depth <= 1:
        print(root)
        for f in files:
            print("   ", f)
print(os.getcwd())
print(os.listdir("."))
