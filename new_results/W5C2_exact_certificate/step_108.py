import subprocess, time
time.sleep(300)
print(subprocess.run(['bash','-c','tail -12 step_03.out'], capture_output=True, text=True).stdout)