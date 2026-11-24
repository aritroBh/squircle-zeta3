import os
d = "../FN2_every_angle_proofs"
print(os.listdir(d))
for f in os.listdir(d):
    p = os.path.join(d, f)
    print(f, os.path.getsize(p))
print("---- my dir ----")
print(os.listdir("."))
for f in ["step_01.py","step_02.py"]:
    print("=====", f, "=====")
    print(open(f).read()[:2000])
