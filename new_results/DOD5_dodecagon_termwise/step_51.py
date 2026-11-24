import os
base = "<project>/research/08_theory/"
D1 = base+"DOD1_dodecagon/"
for f in sorted(os.listdir(D1)):
    if f.startswith("step_") and f.endswith(".py"):
        txt = open(D1+f).read()
        if "Nieint" in txt and "def T3p" in txt:
            print("=====", f)
            i = txt.find("def Nieint")
            print(txt[max(0,i-200):i+2600])
            break
