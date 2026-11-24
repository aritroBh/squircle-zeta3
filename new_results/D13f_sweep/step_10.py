
import os
base = "research/08_theory"
for d in ["D13e_certify"]:
    p = os.path.join(base, d)
    print("==", d, "==")
    for fn in sorted(os.listdir(p)):
        print("  ", fn, os.path.getsize(os.path.join(p,fn)))
