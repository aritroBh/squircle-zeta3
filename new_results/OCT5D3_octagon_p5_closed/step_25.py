base = "<project>/research/08_theory/OCT5D2_octagon_depth_test/"
src = open(base+"verify_oct5d2.py").read()
i = src.find("# ---------- exact rewrite ----------")
j = src.find("# ---------- symbols and controls ----------")
print(src[i:j])
