import os
base = '<project>/research/08_theory'
work = base + '/TOW4_tower_depth8'
# complete_bell region of tow2lib
with open(work + '/tow2lib.py') as f:
    txt = f.read()
i0 = txt.find('def complete_bell')
print(txt[i0:i0+1400])
print("="*80)
with open(work + '/prove_adj56.py') as f:
    print(f.read())
