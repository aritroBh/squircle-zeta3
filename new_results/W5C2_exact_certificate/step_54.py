import subprocess, time
time.sleep(120)
print(subprocess.run(['tail','-3','<project>/research/08_theory/W5C2_exact_certificate/step_01.out'], capture_output=True, text=True).stdout)