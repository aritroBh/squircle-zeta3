code = open('steps/step_07.py').read()
code = code.replace("umax=mp.mpf('5.0')", "umax=mp.mpf('4.6')")
open('steps/step_07.py','w').write(code)
import subprocess, time
t0=time.time()
r=subprocess.run(['<project>/.venv/bin/python','steps/step_07.py','100'],capture_output=True,text=True,timeout=9000)
print("wall: %.1f s"%(time.time()-t0))
print(r.stdout); print("ERR:", r.stderr[-1500:])
