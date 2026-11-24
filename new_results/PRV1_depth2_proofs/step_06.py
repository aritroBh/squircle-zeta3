
import os, subprocess
work = os.getcwd()
print("CWD:", work)
PY='<project>/.venv/bin/python'
r = subprocess.run([PY, 'research/08_theory/CX1_octagon_depth2_proof/verify.py'],
                   capture_output=True, text=True)
print(r.stdout[-2500:])
print('STDERR', r.stderr[-300:])
