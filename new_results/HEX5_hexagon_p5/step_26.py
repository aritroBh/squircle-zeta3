import subprocess, os
print(os.path.exists('result.md'), os.path.getsize('result.md') if os.path.exists('result.md') else 0)
r = subprocess.run(['<project>/.venv/bin/python','step_09.py'],capture_output=True,text=True)
print("STDOUT:", r.stdout); print("STDERR:", r.stderr[-300:])
print(open('result.md').read()[:800])
