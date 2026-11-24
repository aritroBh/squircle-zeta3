import subprocess, time
print(subprocess.run(['tail','-5','<project>/research/08_theory/W5C2_exact_certificate/step_01.out'], capture_output=True, text=True).stdout)
print(subprocess.run(['vm_stat'], capture_output=True, text=True).stdout.splitlines()[:3])