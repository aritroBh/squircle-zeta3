import os
base = '<project>/research/08_theory'
with open(base + '/CX8_tow2_review/report.md') as f:
    txt = f.read()
print("CX8:", len(txt))
print(txt[:7000])
