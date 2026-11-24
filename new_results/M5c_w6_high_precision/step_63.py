import datetime, subprocess, os
PY='<project>/.venv/bin/python'
work='<project>/research/08_theory/M5c_w6_high_precision'
with open(work+'/progress.log','a') as f:
    f.write(f'{datetime.datetime.now().isoformat()} LLL @3300 done (3916s): NO target relation. Best spurious maxcoef 2.3e19-3.0e19, residuals ~1e-3282. Genuine coeffs would exceed ~1e19.5 -> implausible; switch to reserve extras ladder (14,15,13,9,7,6).\n')
# symbol stage for extra pool 14 = (4-3sqrt2)/8  (light, exact)
t0=__import__('time').monotonic()
r=subprocess.run(['nice','-n','10',PY,'-B','w6_solve.py','14'],capture_output=True,text=True,cwd=work,timeout=1700)
print(r.stdout[-1500:]); print('STDERR',r.stderr[-300:])
print('elapsed %.1fs'%(__import__('time').monotonic()-t0))
