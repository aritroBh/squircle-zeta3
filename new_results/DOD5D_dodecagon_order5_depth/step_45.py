import subprocess, time
r = subprocess.run(["<project>/.venv/bin/python", "itab_build.py", "130"],
                   capture_output=True, text=True, timeout=1700)
print(r.stdout)
print(r.stderr[-2000:])
