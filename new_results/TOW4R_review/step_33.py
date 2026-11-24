
import os, glob
base = "<project>"
wts = base + "/research/08_theory/WTS_weight_k_structure/theory"
bs = sorted(glob.glob(wts + "/binteg_*.pkl"), key=lambda p: int(p.split("_")[-1].split(".")[0]))
print("binteg degrees:", [int(p.split("_")[-1].split(".")[0]) for p in bs])
txt = open(base + "/research/08_theory/TOW2_tower_all_depths/fastred.py").read()
i = txt.index("class QSI")
print(txt[i:i+1800])
