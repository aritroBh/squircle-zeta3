import subprocess, time
PY='<project>/research/08_theory/M5c_w6_high_precision'
work='<project>/research/08_theory/M5c_w6_high_precision'
t0=time.monotonic()
r=subprocess.run(['nice','-n','10',PY,'-B','w6c_values_extra.py','14','15'],capture_output=True,text=True,cwd=work,timeout=900)
print(r.stdout[-400:]); print('STDERR',r.stderr[-200:])
print('elapsed %.1fs'%(time.monotonic()-t0))
# also prep symbol stage for {14,15,13}
t0=time.monotonic()
r=subprocess.run(['nice','-n','10',PY,'-B','w6_solve.py','14','15','13'],capture_output=True,text=True,cwd=work,timeout=1700)
print(r.stdout[-350:])
print('solve13 elapsed %.1fs'%(time.monotonic()-t0))
