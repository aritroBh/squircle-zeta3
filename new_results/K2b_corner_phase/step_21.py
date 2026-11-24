
wd = 'research/08_theory/K2b_corner_phase/'
code = open(wd+'saddle2.py').read().replace("return Nsig*(-1)**j*sp.factorial(j-1)","return Nsig*(-1)**(j-1)*sp.factorial(j-1)")
open(wd+'saddle3.py','w').write(code)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python', wd+'saddle3.py'], capture_output=True, text=True)
print(r.stdout); print(r.stderr[-2000:] if r.stderr else '')
