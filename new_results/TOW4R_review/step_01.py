
import os, subprocess, shutil
base = "<project>"
work = os.path.join(base, "research/08_theory/TOW4R_review")
os.makedirs(work, exist_ok=True)
# list the relevant directories
for d in ["research/08_theory", "research"]:
    p = os.path.join(base, d)
    print(d, ":")
    try:
        for x in sorted(os.listdir(p)):
            print("   ", x)
    except Exception as e:
        print("   ERR", e)
