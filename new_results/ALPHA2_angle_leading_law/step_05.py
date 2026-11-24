import subprocess
r = subprocess.run(['cat', '<project>/research/08_theory/ALPHA1_angle_family_across_k/result.md'], capture_output=True, text=True).stdout
print(r[7000:])
