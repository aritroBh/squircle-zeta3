import time, subprocess
time.sleep(420)
print(open('step11b.log').read())
print("running:", subprocess.run(['pgrep','-f','step_11b_d13_lll'],capture_output=True,text=True).stdout.strip() != '')
