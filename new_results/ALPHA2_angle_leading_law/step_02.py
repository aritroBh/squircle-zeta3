import subprocess
print(subprocess.run(['ls', '-la', '<project>/research/08_theory/ALPHA1_angle_family_across_k/'], capture_output=True, text=True).stdout)
