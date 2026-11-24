base = "<project>/research/"
for fn in ["cor6_nf.py", "cor6_ft.py"]:
    with open(base + "08_theory/COR6_true_point_symbols/" + fn) as f:
        print(f.read())
    print("#"*110)
