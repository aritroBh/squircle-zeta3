
import subprocess, os
PY='<project>/.venv/bin/python'
work='research/08_theory/PRV1_depth2_derivative_proof'
os.makedirs(work, exist_ok=True)
os.chdir(work)
r = subprocess.run([PY, 'research/08_theory/CX1_octagon_depth2_proof/verify.py'],
                   capture_output=True, text=True)
print(r.stdout[-3000:])
print('STDERR', r.stderr[-500:])
open('progress.log','w').write("Turn 1: ran CX1 verify.py -> tail below\n"+r.stdout[-1500:]+"\n")
