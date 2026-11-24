base = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
import glob
for f in sorted(glob.glob(base+"step_*.py")):
    n = int(f.split("_")[-1].split(".")[0])
    if n >= 53:
        src = open(f).read()
        if len(src) < 4000:
            print("="*30, f.split("/")[-1], "="*30)
            print(src)
