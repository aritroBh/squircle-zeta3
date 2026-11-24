import os, subprocess, datetime
wd = '<project>/research/08_theory/TOW5R_review'
PY = '<project>/.venv/bin/python'
# delete the copied forced_F910.pkl so depths 9/10 are rebuilt from forced_F78.pkl + arho
os.remove(wd+'/tow5_copy/forced_F910.pkl')
print('forced_F910.pkl deleted from tow5_copy:', not os.path.exists(wd+'/tow5_copy/forced_F910.pkl'))
t0 = datetime.datetime.now()
r = subprocess.run(['nice','-n','10',PY,'forced910.py','9'], cwd=wd+'/tow5_copy', capture_output=True, text=True, timeout=1700)
print(r.stdout)
print('STDERR:', r.stderr[-300:], 'rc:', r.returncode, 'elapsed:', datetime.datetime.now()-t0)
with open(wd+'/progress.log','a') as f:
    f.write(f"{datetime.datetime.now():%H:%M:%S} step_08a: forced910.py 9 (cache deleted, rebuilt) rc={r.returncode} ({(datetime.datetime.now()-t0).seconds}s)\n")
