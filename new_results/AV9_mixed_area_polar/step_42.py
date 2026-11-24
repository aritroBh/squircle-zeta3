
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','pslq2.py','7','300'],capture_output=True,text=True,timeout=1500)
print(r.stdout, r.stderr[-200:])
