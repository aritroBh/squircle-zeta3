import os
base = "<project>"
with open(os.path.join(base, "paper/ieee/README_verification_map.md")) as f:
    rm = f.read()
i = rm.find("$2-L_n=\\sum_{k\\le3}")
if i<0:
    i = rm.find("k\\le3")
print("idx", i)
print(rm[i-500:i+4000])
