import subprocess, os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_02.py','300'],
                   capture_output=True, text=True, timeout=3000)
print(r.stdout[-600:]); print('ERR', r.stderr[-300:])
open('progress.log','a').write('step_02b: extended exact table to N=300, K=6\n')
