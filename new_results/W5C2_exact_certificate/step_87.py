import subprocess
print(subprocess.run(['bash','-c','grep -E "plateau|feed done|DONE|SAVED" step_01b.out'], capture_output=True, text=True).stdout)