import os
base = '<project>/research/08_theory'
with open(base + '/TOW3R_review/result_final.md') as f:
    txt = f.read()
print(len(txt))
print(txt[:9000])
