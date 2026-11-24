import subprocess
print(subprocess.run(['tail','-4','<project>/research/08_theory/ALPHA2_angle_leading_law/kern_old_ext.out'], capture_output=True, text=True).stdout)
