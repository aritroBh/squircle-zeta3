base = "<project>/research/08_theory/"
for f in ["E5c_nonlocal_z5/symeng.py", "OCT5D2_octagon_depth_test/oct5d2_engine.py"]:
    print("="*100); print(f); print("="*100)
    with open(base+f) as fh:
        print(fh.read())
