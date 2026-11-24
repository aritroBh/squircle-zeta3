base = "<project>/research/08_theory/"
for f in ["gengine.py", "oct5d2_engine.py", "symengine.py"]:
    print("="*100)
    print("FILE:", f)
    with open(base + "OCT5D2_octagon_depth_test/" + f) as fh:
        print(fh.read())
