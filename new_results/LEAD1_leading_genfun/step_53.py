import subprocess, os, time
os.chdir('<project>/research/08_theory/LEAD1_leading_genfun')
PY = '<project>/.venv/bin/python'
t0 = time.time()
r = subprocess.run(['nice','-n','10',PY,'-B','verify_lead1.py'], capture_output=True, text=True, timeout=1500)
print(r.stdout)
print('STDERR tail:', r.stderr[-500:])
print('rc=', r.returncode, f'{time.time()-t0:.0f}s')
