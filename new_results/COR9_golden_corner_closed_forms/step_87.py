import mpmath as mp
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
mp.mp.dps = 100
for tag in ("b2sq", "gol1", "gol3"):
    v = mp.mpf(open(WD+f"{tag}_M3_hpl_1150.txt").read())
    print(tag, "M3 =", mp.nstr(v, 60))
