import os
base = '<project>/research/08_theory'
tow3 = base + '/TOW3_tower_depth6'
for sub in ['logic56', 'rerun', 'w5lens']:
    print(sub, ':', sorted(os.listdir(tow3 + '/claude_screen/' + sub)))
print("="*80)
t3r = base + '/TOW3R_review'
print(sorted(os.listdir(t3r)))
