import subprocess, time, os
PY='<project>/.venv/bin/python'
work='<project>/research/08_theory/M5c_w6_high_precision'
r=subprocess.run(['nice','-n','10',PY,'-B',work+'/w6c_values_extra.py','14','15','13','9','7'],capture_output=True,text=True,cwd=work,timeout=900)
print(r.stdout[-200:]); print('STDERR',r.stderr[-150:])
t0=time.monotonic()
r=subprocess.run(['nice','-n','10',PY,'-B',work+'/w6_solve.py','14','15','13','9','7','6'],capture_output=True,text=True,cwd=work,timeout=1700)
print(r.stdout[-300:]); print('STDERR',r.stderr[-150:])
print('solve6 elapsed %.1fs'%(time.monotonic()-t0))
