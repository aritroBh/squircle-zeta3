import subprocess, datetime
wd = '<project>/research/08_theory/TOW5R_review'
PY = '<project>/.venv/bin/python'
cases = [(9,0),(9,1),(9,2),(9,3),(9,4),(10,0),(10,1),(10,2),(10,3),(10,4),(10,5)]
t00 = datetime.datetime.now()
allok = True
for (dl, s) in cases:
    t0 = datetime.datetime.now()
    r = subprocess.run(['nice','-n','10',PY,'prove_adj910_fast.py',str(dl),str(s)], cwd=wd+'/tow5_copy', capture_output=True, text=True, timeout=1700)
    ok = (r.returncode == 0) and ('numerator == 0 : True' in r.stdout)
    allok &= ok
    print(f'({dl},{s}): rc={r.returncode} ok={ok} ({(datetime.datetime.now()-t0).seconds}s) {r.stdout.strip().splitlines()[-1] if r.stdout else ""}', flush=True)
with open(wd+'/progress.log','a') as f:
    f.write(f"{datetime.datetime.now():%H:%M:%S} step_10: prove_adj910_fast.py all 11 identities rc0+True: {allok} ({(datetime.datetime.now()-t00).seconds}s total)\n")
print('ALL 11:', allok)
