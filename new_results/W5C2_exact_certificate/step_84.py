import subprocess
out = subprocess.run(['bash','-c','grep -E "ppb|step_01b pp" step_01b.out | tail -25'], capture_output=True, text=True).stdout
print(out)