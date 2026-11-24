base = "<project>/research/08_theory/"
for f in ["OCT5D_octagon_order5_depth/polygon_engine.py", "OCT5D_octagon_order5_depth/oct_engine.py"]:
    print("="*100); print(f); print("="*100)
    with open(base+f) as fh:
        print(fh.read())
