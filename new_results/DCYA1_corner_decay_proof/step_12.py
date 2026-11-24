import subprocess
d='<project>/research/08_theory/D10c_corner_general/'
print("=== files ===")
print(subprocess.run(['ls','-la',d], capture_output=True, text=True).stdout)
# Find the main result/readme file
for f in ['result.md','proof.md','README.md','BRIEF.md']:
    import os
    if os.path.exists(d+f):
        print(f"\n=== {f} ===")
        print(subprocess.run(['cat', d+f], capture_output=True, text=True).stdout[:6000])
