import os
base = '<project>/research/08_theory'
with open(base + '/TOW2R_review/result.md') as f:
    txt = f.read()
print("TOW2R:", len(txt))
print(txt[:5500])
