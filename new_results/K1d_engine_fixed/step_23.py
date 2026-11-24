s = open('step_43.py').read()
s = s.replace("rel = mp.pslq(basis_num, tol=mp.mpf('1e-100'), maxcoeff=mp.mpf('1e12'), maxsteps=10000)",
              "rel = mp.pslq(basis_num, tol=mp.mpf('1e-95'), maxcoeff=10**12, maxsteps=10000)")
open('step_43.py','w').write(s)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_43.py'],
                   capture_output=True, text=True, timeout=900)
print(r.stdout)
print("STDERR:", r.stderr[-3000:] if r.stderr else "")
