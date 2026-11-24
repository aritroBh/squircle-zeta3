import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_40.py'],
                   capture_output=True, text=True, timeout=600)
print(r.stdout)
print(r.stderr)
with open('progress.log','a') as f:
    f.write("step_40: rebuilt E3_fixed.pkl with 1/4 fix. "+r.stdout.strip().replace("\n"," | ")+"\n")
