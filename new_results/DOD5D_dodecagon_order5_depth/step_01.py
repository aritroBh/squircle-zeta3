import os, glob
base = "<project>/research/08_theory/"
for d in ["OCT5D_octagon_order5_depth", "OCT5D2_octagon_depth_test", "M4c_w5_closed_form", "DOD5_dodecagon_termwise", "DOD1_dodecagon"]:
    p = os.path.join(base, d)
    print("="*80)
    print(d, "exists:", os.path.isdir(p))
    for f in sorted(glob.glob(p+"/**", recursive=True))[:40]:
        print("  ", f.replace(p, ""), os.path.getsize(f) if os.path.isfile(f) else "<dir>")
