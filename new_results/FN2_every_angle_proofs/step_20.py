import os
base = "<project>"
ah3 = os.path.join(base, "research/08_theory/AH3_fn_any_angle")
with open(os.path.join(ah3, "step_29.py")) as f:
    print(f.read()[:6000])
