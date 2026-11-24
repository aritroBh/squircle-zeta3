import subprocess, os, time, datetime
PY='<project>/.venv/bin/python'
work='<project>/research/08_theory/M5c_w6_high_precision'
t0=time.monotonic()
r=subprocess.run(['nice','-n','10',PY,'-B','w6c_values_extra.py','14'],capture_output=True,text=True,cwd=work,timeout=1700)
print(r.stdout[-600:]); print('STDERR',r.stderr[-200:])
print('elapsed %.1fs'%(time.monotonic()-t0))
with open(work+'/progress.log','a') as f:
    f.write(f'{datetime.datetime.now().isoformat()} extra pool14=(4-3sqrt2)/8: symbol stage rank60 kernel90 target-in-span; Li6 leaf+atom values at 3560 digits saved.\n')
