import time
time.sleep(240)
print(open('step09e.log').read())
print("running:", end=" ")
import subprocess
print(subprocess.run(['pgrep','-f','step_09e_d13_kernel'],capture_output=True,text=True).stdout.strip() != '')
