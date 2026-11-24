import os, glob
base = "<project>/research/08_theory/"
print(os.getcwd())
print(sorted(os.listdir(base)))
for d in sorted(os.listdir(base)):
    p = os.path.join(base, d)
    if os.path.isdir(p) and d.startswith("DOD"):
        print("\n==", d)
        for f in sorted(os.listdir(p)):
            print("   ", f)
