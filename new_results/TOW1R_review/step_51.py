import os
work = "<project>/research/08_theory/TOW1R_review/"
print(open(work+"progress.log").read())
print("--- files ---")
for f in sorted(os.listdir(work)):
    if f.startswith("step_") or f in ("myengine.py","result.md"):
        print("  ",f, os.path.getsize(work+f))
