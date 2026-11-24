import os
base = '<project>/research/08_theory'
tow3 = base + '/TOW3_tower_depth6'
with open(tow3 + '/tow2lib.py') as f:
    txt = f.read()
print(txt[6000:9056])
print("="*80)
print(sorted(os.listdir(tow3 + '/claude_screen')))
