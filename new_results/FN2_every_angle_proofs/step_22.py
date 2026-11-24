import os
wd = "<project>/research/08_theory/FN2_every_angle_proofs"
os.chdir(wd)
# Check what these files are
for fn in ["step_01.py","step_02.py","step_03.py"]:
    with open(fn) as f:
        head = f.read()[:400]
    print("="*20, fn, "="*20)
    print(head)
import subprocess
print(subprocess.run(["ls","-la","research"],capture_output=True,text=True).stdout)
