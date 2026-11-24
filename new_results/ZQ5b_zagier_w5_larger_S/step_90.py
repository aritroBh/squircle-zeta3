import time, subprocess
time.sleep(420)
print(open('step12.log').read())
print("running:", subprocess.run(['pgrep','-f','step_12_d13_red_lll'],capture_output=True,text=True).stdout.strip() != '')
