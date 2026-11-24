
import os
base = "<project>"
for d in ["TOW2R_review", "TOW3R_review"]:
    p = os.path.join(base, "research/08_theory", d)
    print("="*20, d)
    for x in sorted(os.listdir(p)):
        print(f"  {os.path.getsize(os.path.join(p,x)):>9}  {x}")
import glob
for pat in ["research/08_theory/CX8*", "research/08_theory/CX9*"]:
    for p in glob.glob(os.path.join(base, pat)):
        print(p)
        if os.path.isdir(p):
            for x in sorted(os.listdir(p))[:40]:
                print(f"   {os.path.getsize(os.path.join(p,x)):>9}  {x}")
