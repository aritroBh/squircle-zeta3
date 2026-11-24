import subprocess, time, os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
time.sleep(60)
print(subprocess.run(['tail','-3','step_01.out'], capture_output=True, text=True).stdout)
print(subprocess.run(['bash','-c','ps -A -o pid,time,command 2>&1 | head -6'], capture_output=True, text=True).stdout)