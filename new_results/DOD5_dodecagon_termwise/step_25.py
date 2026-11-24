import os
base = "<project>/research/08_theory/"
D2 = base+"DOD2_dodecagon_finish/"
for f in sorted(os.listdir(D2)):
    if f.startswith("step_") and f.endswith(".py"):
        txt = open(D2+f).read()
        if "Xi_verified" in txt and "dump" in txt:
            print("=====", f)
            print(txt[:4500])
