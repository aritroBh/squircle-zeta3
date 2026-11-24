import subprocess
r = subprocess.run(["<project>/.venv/bin/python", "-u", "dod5d_assembly.py", "itab"],
                   capture_output=True, text=True, timeout=600)
print(r.stdout); print(r.stderr[-300:])
