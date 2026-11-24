import subprocess, time
time.sleep(500)
print(subprocess.run(['tail','-8','<project>/research/08_theory/W5C2_exact_certificate/step_01.out'], capture_output=True, text=True).stdout)