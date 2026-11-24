import subprocess, time
time.sleep(360)
print(subprocess.run(['bash','-c','cat step_03_test.out'], capture_output=True, text=True).stdout)