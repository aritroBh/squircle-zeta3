import os
base = "<project>/research/08_theory/"
D1 = base+"DOD1_dodecagon/"
for f in sorted(os.listdir(D1)):
    if f.startswith("step_") and f.endswith(".py"):
        txt = open(D1+f).read()
        if "atoms620" in txt and "dump" in txt:
            print("=====", f)
            print(txt[:5000])
