import subprocess, time
time.sleep(180)
print(subprocess.run(['bash','-c','tail -8 step_02.out'], capture_output=True, text=True).stdout)