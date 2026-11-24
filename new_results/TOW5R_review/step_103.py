import subprocess, datetime
wd = '<project>/research/08_theory/TOW5R_review'
PY = '<project>/.venv/bin/python'
t0 = datetime.datetime.now()
r = subprocess.run(['nice','-n','10',PY,'verify_tow5r.py'], cwd=wd, capture_output=True, text=True, timeout=1700)
print(r.stdout)
print('rc:', r.returncode, 'STDERR tail:', r.stderr[-300:], 'elapsed:', datetime.datetime.now()-t0)
open(wd+'/verify_tow5r.log','w').write(r.stdout + '\nSTDERR:\n' + r.stderr)
with open(wd+'/progress.log','a') as f:
    f.write(f"{datetime.datetime.now():%H:%M:%S} step_23: verify_tow5r.py full run rc={r.returncode} ({(datetime.datetime.now()-t0).seconds}s)\n")
