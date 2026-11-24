import subprocess
for f in ['step_03.py', 'step_04.py']:
    r = subprocess.run(['python', f], capture_output=True, text=True, timeout=1200)
    print("="*20, f, "rc=", r.returncode)
    print(r.stdout[-3000:])
    if r.returncode != 0: print(r.stderr[-1500:])
