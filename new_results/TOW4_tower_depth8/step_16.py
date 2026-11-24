import os
base = '<project>/research/08_theory'
with open(base + '/CX8_tow2_review/report.md') as f:
    txt = f.read()
print(txt[11300:18885])
