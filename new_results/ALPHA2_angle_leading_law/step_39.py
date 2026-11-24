import subprocess
print(subprocess.run(['tail','-3','<project>/research/08_theory/ALPHA2_angle_leading_law/kern_new.out'], capture_output=True, text=True).stdout)
