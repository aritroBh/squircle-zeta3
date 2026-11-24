base = "<project>/research/08_theory/"
for f in ["OCT5D_octagon_order5_depth/result.md"]:
    with open(base+f) as fh:
        print(fh.read())
