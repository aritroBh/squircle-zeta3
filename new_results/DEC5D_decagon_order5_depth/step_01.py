
import os
base = "<project>/research/08_theory/"
for d in ["OCT5D_octagon_order5_depth", "OCT5D2_octagon_depth_test", "M4c_w5_closed_form"]:
    p = os.path.join(base, d)
    print("="*80)
    print(d, os.path.isdir(p))
    if os.path.isdir(p):
        for f in sorted(os.listdir(p)):
            print("   ", f, os.path.getsize(os.path.join(p,f)))
print("="*80)
print("cwd:", os.getcwd())
print(os.listdir("."))
