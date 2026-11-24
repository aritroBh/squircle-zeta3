import os, subprocess
base = '<project>/research/08_theory'
tow3 = base + '/TOW3_tower_depth6'
for s in ['step_63.py','step_65.py','step_68.py']:
    out = subprocess.run(['grep','-l','kern_5\|cubic\|40/40\|40 / 40', tow3+'/'+s], capture_output=True, text=True)
    print(s, out.stdout.strip())
print("="*60)
print(open(tow3 + '/step_63.py').read()[:5000])
