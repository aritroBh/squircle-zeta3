import os

base = "<project>/research/08_theory/"
for d in ["OCT5D2_octagon_depth_test", "OCT5D_octagon_order5_depth", "M4e_w5_lchi", "M4f_w5_extras"]:
    p = os.path.join(base, d)
    print("="*80)
    print(d)
    for root, dirs, files in os.walk(p):
        for f in files:
            fp = os.path.join(root, f)
            print(f"  {fp.replace(base,'')}  ({os.path.getsize(fp)} bytes)")
