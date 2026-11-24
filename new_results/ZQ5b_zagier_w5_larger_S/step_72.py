import time, subprocess
time.sleep(300)
print(open('step09e.log').read())
print("running:", subprocess.run(['pgrep','-f','step_09e_d13_kernel'],capture_output=True,text=True).stdout.strip() != '')
