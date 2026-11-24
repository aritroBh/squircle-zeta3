wd = '<project>/research/08_theory/TOW5R_review'
s = open(wd+'/step_06.py').read()
s = s.replace("from mpmath import mp, mpc, power\nmp.mp.dps = 40", "from mpmath import mp, mpc, power\nmp.dps = 40")
open(wd+'/step_06.py','w').write(s)
import subprocess, datetime
PY = '<project>/.venv/bin/python'
t0 = datetime.datetime.now()
r = subprocess.run(['nice','-n','10',PY,'step_06.py'], cwd=wd, capture_output=True, text=True, timeout=1700)
print(r.stdout)
print('STDERR tail:', r.stderr[-400:], 'rc:', r.returncode, 'elapsed:', datetime.datetime.now()-t0)
with open(wd+'/progress.log','a') as f:
    f.write(f"{datetime.datetime.now():%H:%M:%S} step_06 rerun: rc={r.returncode} ({(datetime.datetime.now()-t0).seconds}s)\n")
