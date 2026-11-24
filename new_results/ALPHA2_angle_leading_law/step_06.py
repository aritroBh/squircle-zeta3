import subprocess
for f in ['step_01.py','step_02.py','step_03.py','step_04.py']:
    print('='*20, f)
    print(subprocess.run(['cat', f'<project>/research/08_theory/ALPHA2_angle_leading_law/{f}'], capture_output=True, text=True).stdout)
