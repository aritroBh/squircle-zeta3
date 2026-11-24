import subprocess, time
time.sleep(240)
print(subprocess.run(['bash','-c','tail -8 step_05.out'], capture_output=True, text=True).stdout)