
import os
base = "<project>"
d = os.path.join(base,"research/08_theory/GRW1_diamond_growth_borel")
for f in sorted(os.listdir(d)):
    if f.endswith(".py"):
        src = open(os.path.join(d,f)).read()
        print("="*30, f, " len:", len(src))
        # print first 40 lines of each
        print("\n".join(src.splitlines()[:25]))
