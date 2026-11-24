import subprocess, time, os
PY='<project>/.venv/bin/python'
work='<project>/research/08_theory/M5c_w6_high_precision'
time.sleep(20)
t0=time.monotonic()
r=subprocess.run(['nice','-n','10',PY,'-B',work+'/w6c_values_extra.py','14','15'],capture_output=True,text=True,cwd=work,timeout=900)
print(r.stdout[-400:]); print('STDERR',r.stderr[-300:])
print('elapsed %.1fs'%(time.monotonic()-t0))
print(os.path.exists(work+'/w6_atomvals_3500_x14_15.pkl'), os.path.exists(work+'/w6_leaves_extra_w6x14_15.pkl'))
