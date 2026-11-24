import subprocess
r = subprocess.run(['<project>/.venv/bin/python','verify_pwr1.py'],
                   capture_output=True, text=True, timeout=280, cwd='.')
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "")
