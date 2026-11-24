import subprocess, pathlib
PY = '<project>/.venv/bin/python'
print(subprocess.run(['bash','-c','grep -n "2\\*\\*(d-w)\\|2\\*\\*(w-d)\\|cS = P.pop" rul1.py'], capture_output=True, text=True).stdout)
