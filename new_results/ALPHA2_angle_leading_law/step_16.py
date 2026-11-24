import subprocess
A1='<project>/research/08_theory/ALPHA1_angle_family_across_k'
for f in ['step_08.py','step_31.py','step_32.py']:
    print('='*25, f)
    print(subprocess.run(['cat', f'{A1}/{f}'], capture_output=True, text=True).stdout)
