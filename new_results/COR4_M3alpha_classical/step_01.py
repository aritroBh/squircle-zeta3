import os
base = "<project>/research/"
for root, dirs, files in os.walk(base):
    if any(k in root for k in ["G1_zeta3_any_angle", "M4c_w5_closed_form", "COR4", "planner"]):
        print(root)
        for f in files[:50]:
            print("   ", f)
