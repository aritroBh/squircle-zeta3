import os
base = "<project>"
with open(os.path.join(base, "research/08_theory/AH3_fn_any_angle/result.md")) as f:
    ah3 = f.read()
print("AH3 result.md length:", len(ah3))
print(ah3[:12000])
