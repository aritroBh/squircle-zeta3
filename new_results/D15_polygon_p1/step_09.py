import os
for d in ["D10c_corner_general", "D1b_corner_identify", "D1_corner_level", "D10b_octagon"]:
    p = f"<project>/research/08_theory/{d}"
    if os.path.isdir(p):
        print(d, os.listdir(p)[:30])
    else:
        print(d, "NOT A DIR")
