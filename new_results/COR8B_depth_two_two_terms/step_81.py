
content = open('steps/step_11.py').read()
content = content.replace("mp.mpf(Fr(t[1]))", "mpf_f(Fr(t[1]))")
content = content.replace("mp.mpf(cf)", "mpf_f(cf)")
content = content.replace("import sympy as sp", "import sympy as sp\ndef mpf_f(fr): return mp.mpf(fr.numerator)/mp.mpf(fr.denominator)")
open('steps/step_11.py','w').write(content)
import subprocess, time
t=time.time()
r = subprocess.run(['nice','-n','15','<project>/.venv/bin/python','steps/step_11.py'],capture_output=True,text=True,timeout=1400)
print(r.stdout[-4000:]); print(r.stderr[-3000:]); print("elapsed", time.time()-t)
