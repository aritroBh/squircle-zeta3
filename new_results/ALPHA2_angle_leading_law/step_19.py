import subprocess
print(subprocess.run(['cat', '<project>/research/08_theory/ALPHA2_angle_leading_law/cwred.py'], capture_output=True, text=True).stdout)
