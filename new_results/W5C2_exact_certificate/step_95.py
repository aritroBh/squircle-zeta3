import subprocess, time
time.sleep(400)
print(subprocess.run(['bash','-c','grep -E "ALL CHECKS|DONE|step_02 (pm|mp|mm)" step_02.out | tail -12'], capture_output=True, text=True).stdout)