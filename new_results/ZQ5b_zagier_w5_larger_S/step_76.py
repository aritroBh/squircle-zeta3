import time, subprocess
time.sleep(240)
print(open('step10.log').read())
print("running:", subprocess.run(['pgrep','-f','step_10_d13_vals'],capture_output=True,text=True).stdout.strip() != '')
