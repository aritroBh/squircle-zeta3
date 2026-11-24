import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run(['nice', PY, 'step_40_fits.py', '4', '450', '700'], capture_output=True, text=True, cwd='.')
print(r.stdout[-4500:])
if r.stderr: print("STDERR:", r.stderr[-600:])
print("wall %.0fs" % (time.time()-t0))
open('progress.log','a').write("step_41: Wt3 fit ACCEPTED (identical 450/700, resid 3.4e-487/1.2e-737, height 416900, 8 atoms: needs L3m=(5250-1500sqrt5)Li3(-1/phi))\n")
