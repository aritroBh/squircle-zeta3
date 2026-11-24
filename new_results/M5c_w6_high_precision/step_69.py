import subprocess, time, os
PY='<project>/.venv/bin/python'
work='<project>/research/08_theory/M5c_w6_high_precision'
t0=time.monotonic()
r=subprocess.run(['nice','-n','10',PY,'-B',work+'/w6_solve.py','14','15','13'],capture_output=True,text=True,cwd=work,timeout=1700)
print(r.stdout[-350:]); print('STDERR',r.stderr[-200:])
print('elapsed %.1fs'%(time.monotonic()-t0))
