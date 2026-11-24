import subprocess, os, time, datetime, signal
PY='<project>/.venv/bin/python'
work='<project>/research/08_theory/M5c_w6_high_precision'
# prep {14,15,9} symbol stage + values (light)
t0=time.monotonic()
r=subprocess.run(['nice','-n','10',PY,'-B',work+'/w6_solve.py','14','15','9'],capture_output=True,text=True,cwd=work,timeout=1700)
print(r.stdout[-280:])
r=subprocess.run(['nice','-n','10',PY,'-B',work+'/w6c_values_extra.py','14','15','9'],capture_output=True,text=True,cwd=work,timeout=900)
print(r.stdout[-200:])
print('prep elapsed %.1fs'%(time.monotonic()-t0))
# kill the {14,15,13} run
pid=int(open(work+'/w6c_lll_x14_15_13_3300.pid').read())
os.kill(pid, signal.SIGTERM)
time.sleep(3)
try: os.kill(pid,0); print('still alive')
except OSError: print('killed', pid)
with open(work+'/progress.log','a') as f:
    f.write(f'{datetime.datetime.now().isoformat()} Duplication identities verified @3500 digits: Li6(b^2)=32(Li6(b)+Li6(-b)) (8.5e-3610), Li6(b^4)=1024 base combo (2.8e-3608), Li6(-1/r2)=Li6(1/2)/32-Li6(1/r2) (5.1e-3611) -> extras 13,7,6 are exact base combos, SKIPPED. Killed w6x14_15_13 run (span == w6x14_15 span). Next: w6x14_15_9.\n')
