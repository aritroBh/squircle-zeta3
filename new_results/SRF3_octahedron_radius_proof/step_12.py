
import os
base = "<project>"
for f in ["research/08_theory/GRW1_diamond_growth_borel",
          "research/08_theory/R2_rhombus_ab",
          "research/08_theory/D16_pure_weight"]:
    p = os.path.join(base, f)
    print("==", f)
    if os.path.isdir(p):
        for x in sorted(os.listdir(p)):
            print("   ", x)
    print()
