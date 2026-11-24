base = "<project>/research/08_theory/"
for f in ["OCT5D2_octagon_depth_test/gengine.py", "OCT5D2_octagon_depth_test/symengine.py"]:
    print("="*100); print(f); print("="*100)
    with open(base+f) as fh:
        print(fh.read())
