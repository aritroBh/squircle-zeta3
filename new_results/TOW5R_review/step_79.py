import subprocess, datetime
wd = '<project>/research/08_theory/TOW5R_review'
PY = '<project>/.venv/bin/python'
t00 = datetime.datetime.now()
allok = True
for dl in [9, 10]:
    t0 = datetime.datetime.now()
    r = subprocess.run(['nice','-n','10',PY,'prove_E3b910.py',str(dl)], cwd=wd+'/tow5_copy', capture_output=True, text=True, timeout=1700)
    ok = (r.returncode == 0) and ('numerator == 0 : True' in r.stdout)
    allok &= ok
    print(f'E3({dl}): rc={r.returncode} ok={ok} ({(datetime.datetime.now()-t0).seconds}s) {r.stdout.strip().splitlines()[-1] if r.stdout else ""}', flush=True)
with open(wd+'/progress.log','a') as f:
    f.write(f"{datetime.datetime.now():%H:%M:%S} step_12: prove_E3b910.py 9,10 rc0+True: {allok} ({(datetime.datetime.now()-t00).seconds}s)\n")
print('E3 BOTH:', allok)
