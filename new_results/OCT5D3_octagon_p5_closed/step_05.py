import os
print("CWD:", os.getcwd())
print(os.listdir("."))
base = "<project>/research/08_theory/"
with open(base + "OCT5D_octagon_order5_depth/result.md") as f:
    print(f.read())
