import os
print("CWD:", os.getcwd())
base = "<project>/research/08_theory/"
for root, dirs, files in os.walk(os.path.join(base, "K4e_tetrakis_d3")):
    print(root)
    for f in files:
        p = os.path.join(root, f)
        print("   ", f, os.path.getsize(p))
print()
print(os.listdir(os.path.join(base, "K4g_tetrakis_z3")))
print()
print("cwd listing:", os.listdir("."))
