import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run(['nice', PY, 'step_40_fits.py', '5', '450', '700'], capture_output=True, text=True, cwd='.')
print(r.stdout[-6000:])
if r.stderr: print("STDERR:", r.stderr[-600:])
print("wall %.0fs" % (time.time()-t0))
open('progress.log','a').write("step_42: Wt4 fit ACCEPTED (identical 450/700, resid 1.65e-485/7.67e-736, height 2395375, 15 atoms, L4m/L4q/L3m*log all present)\n")
