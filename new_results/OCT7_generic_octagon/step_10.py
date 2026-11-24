import os, re
base="research/08_theory"
o5=os.path.join(base,"OCT5_octagon_exact")
# grep for E3 construction in OCT5 steps
for f in sorted(os.listdir(o5)):
    if f.startswith("step_") and f.endswith(".py"):
        src=open(os.path.join(o5,f)).read()
        if "E3" in src and ("monom" in src or "edge" in src or "A0" in src or "kappa" in src or "kap" in src):
            print("###", f, "len", len(src))
