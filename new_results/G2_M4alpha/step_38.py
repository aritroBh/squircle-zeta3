import os, re
base = "<project>/research/"
t = open(base + "08_theory/M4e_w5_lchi/result.md").read()
i = t.find("0.0060279931")
print(t[max(0,i-400):i+300])
for d in os.listdir(base + "08_theory"):
    if d.startswith("M4a") and os.path.isdir(base + "08_theory/" + d):
        print("M4a dir:", d, sorted(os.listdir(base + "08_theory/" + d))[:40])
