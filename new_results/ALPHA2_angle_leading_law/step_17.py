import subprocess
for f in ['cwfield.py','cwbuild.py']:
    print('='*25, f)
    print(subprocess.run(['cat', f'<project>/research/08_theory/ALPHA2_angle_leading_law/{f}'], capture_output=True, text=True).stdout)
