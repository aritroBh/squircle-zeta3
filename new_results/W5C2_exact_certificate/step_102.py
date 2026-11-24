import subprocess, time
time.sleep(430)
print(subprocess.run(['bash','-c','tail -6 step_03.out'], capture_output=True, text=True).stdout)