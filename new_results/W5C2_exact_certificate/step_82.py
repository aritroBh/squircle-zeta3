import subprocess, time
time.sleep(420)
print(subprocess.run(['tail','-10','<project>/research/08_theory/W5C2_exact_certificate/step_01b.out'], capture_output=True, text=True).stdout)