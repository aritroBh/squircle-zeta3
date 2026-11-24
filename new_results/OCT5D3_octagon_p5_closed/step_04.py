base = "<project>/research/08_theory/"
with open(base + "OCT5D2_octagon_depth_test/gengine.py") as fh:
    src = fh.read()
print(src[-2000:])
print("="*100)
with open(base + "OCT5D2_octagon_depth_test/verify_oct5d2.py") as fh:
    print(fh.read())
