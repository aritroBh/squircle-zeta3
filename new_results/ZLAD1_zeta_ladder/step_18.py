import subprocess, os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
# quick sanity at small N first
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_02.py','8'],
                   capture_output=True, text=True, timeout=1200)
print(r.stdout[-3000:])
print(r.stderr[-1500:])
