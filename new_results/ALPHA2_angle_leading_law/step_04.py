import subprocess, os
os.makedirs('<project>/research/08_theory/ALPHA2_angle_leading_law', exist_ok=True)
print(subprocess.run(['ls', '<project>/research/08_theory/ALPHA2_angle_leading_law'], capture_output=True, text=True).stdout)
