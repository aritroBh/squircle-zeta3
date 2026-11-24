import subprocess, os
d='research/08_theory/D13h_final_cell'
r = subprocess.run(['<home>)/.venv/bin/python','step_10_J192.py'],
                   cwd=d, capture_output=True, text=True)
print(r.stdout[-3200:])
print("STDERR:", r.stderr[-500:] if r.stderr else "(none)")
