import os, shutil, subprocess, sys
BASE = '<project>'
OUT = os.path.join(BASE,'paper/squircle_note/outputs')
for f in ['m4_tails_exact.pkl','m3_tails_exact.pkl']:
    shutil.copy(os.path.join(OUT,f), f)
r = subprocess.run([sys.executable, 'step_quad_pieces.py'], capture_output=True, text=True, timeout=2400)
print(r.stdout)
print('ERR', r.stderr[-1000:])
import datetime
with open('progress.log','a') as fh:
    fh.write(f"[{datetime.datetime.now().isoformat()}] step_quad_pieces.py: W5,W4'',W3^4,T,M4 quadrature at 80 dps; sum-M4 residual printed.\n")
