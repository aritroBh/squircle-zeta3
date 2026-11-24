import os
base = "<project>"
ah3 = os.path.join(base, "research/08_theory/AH3_fn_any_angle")
for fn in ["step_15.py"]:
    with open(os.path.join(ah3, fn)) as f:
        print("="*30, fn, "="*30)
        print(f.read()[:5000])
