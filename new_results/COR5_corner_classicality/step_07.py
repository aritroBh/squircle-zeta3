import os
base = "<project>/research"
cor4 = os.path.join(base, "08_theory/COR4_M3alpha_classical")
with open(os.path.join(cor4, "cor4_sym.py")) as f:
    txt = f.read()
# find truncation point: print chars 16000 onwards? LEN?
print("LEN:", len(txt))
# print the middle section that was truncated
print(txt[11500:17000])
