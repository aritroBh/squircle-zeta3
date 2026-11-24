s = open('step_45.py').read()
s = s.replace('print("PhiA terms:", len(sp.Poly(PhiA, Jx).terms()) if PhiA!=0 else 0)\n', '')
open('step_45.py','w').write(s)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_45.py'],
                   capture_output=True, text=True, timeout=1800)
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "")
