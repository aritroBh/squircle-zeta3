import subprocess, os
PY = '<project>/.venv/bin/python'
a04dir = '<project>/research/10_frontier/W1_corner_family/A_function_level/'
# run a04 in its own dir (it has the numeric control built in)
r = subprocess.run(['nice','-n','10',PY,'a04_rfamily.py'], capture_output=True, text=True, cwd=a04dir, timeout=1500)
print(r.stdout[-3000:])
print(r.stderr[-1000:] if r.returncode else "RC=0")
