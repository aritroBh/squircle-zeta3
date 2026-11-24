import subprocess, datetime
wd = '<project>/research/08_theory/TOW5R_review'
PY = '<project>/.venv/bin/python'
t0 = datetime.datetime.now()
r = subprocess.run(['nice','-n','10',PY,'verify_L.py'], cwd=wd+'/tow5_copy', capture_output=True, text=True, timeout=1700)
print(r.stdout)
print('STDERR tail:', r.stderr[-300:], 'rc:', r.returncode, 'elapsed:', datetime.datetime.now()-t0)
with open(wd+'/progress.log','a') as f:
    f.write(f"{datetime.datetime.now():%H:%M:%S} step_07a: verify_L.py rerun (L cache == tow2lib L(r), r=1..10, symbolic) rc={r.returncode} ({(datetime.datetime.now()-t0).seconds}s)\n")
