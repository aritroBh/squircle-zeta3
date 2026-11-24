import subprocess, time
time.sleep(660)
print(subprocess.run(['bash','-c','tail -8 step_03.out'], capture_output=True, text=True).stdout)