
import os
base = "research/08_theory/D13e_certify"
for fn in ["dbg1.py", "step_30.py"]:
    print("="*20, fn, "="*20)
    with open(os.path.join(base, fn)) as f:
        print(f.read())
