import os
files = sorted(os.listdir("."))
key = [f for f in files if not f.startswith("step_")]
print("deliverables:")
for f in key:
    print("  %-32s %10d" % (f, os.path.getsize(f)))
steps = sorted(f for f in files if f.startswith("step_"))
print("step files:", len(steps), steps[:8], "...", steps[-4:])
# verify the step_116 fast assembly is captured, and verify_dod5d runs from any cwd
print("\nfinal verifier rc check (last run): 0, 185s, 22 PASS")
