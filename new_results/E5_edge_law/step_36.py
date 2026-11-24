import subprocess
r=subprocess.run(['<project>/.venv/bin/python','steps/step_06.py'],capture_output=True,text=True,timeout=900)
print("OUT:", r.stdout[:2000])
print("ERR FULL:", r.stderr)
