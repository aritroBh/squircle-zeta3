import subprocess, sys
r = subprocess.run(['nice', sys.executable, 'step_03_watson.py'],capture_output=True,text=True,timeout=600,cwd='../DCY2_decay_first_correction')
print(r.stdout); print(r.stderr[-200:] if r.stderr else "")
open('progress.log','a').write("step: verify_dcy2r.py PASS 132s; step_03_watson C_rho rerun\n")
