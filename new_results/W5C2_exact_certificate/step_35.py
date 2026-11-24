import subprocess, time
print(subprocess.run(['tail','-4','step_01.out'], capture_output=True, text=True).stdout)