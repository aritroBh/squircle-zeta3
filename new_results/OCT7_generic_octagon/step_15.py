import os
base="research/08_theory"
o6=os.path.join(base,"OCT6_octagon_closed")
for f in sorted(os.listdir(o6)):
    if f.startswith("step_") and f.endswith(".py"):
        src=open(os.path.join(o6,f)).read()
        if "mixed_engine" in src or "INT(" in src or "MIX(" in src:
            print("###", f, "len", len(src))
