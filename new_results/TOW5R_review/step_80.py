import subprocess, datetime, os, time
wd = '<project>/research/08_theory/TOW5R_review'
PY = '<project>/.venv/bin/python'
t0 = datetime.datetime.now()
r = subprocess.run(['nice','-n','10',PY,'verify_W7.py'], cwd=wd+'/tow5_copy', capture_output=True, text=True, timeout=1700)
print(r.stdout)
print('rc:', r.returncode, 'stderr:', r.stderr[-200:], 'elapsed:', datetime.datetime.now()-t0)
# provenance of kern files
wts = '<project>/research/08_theory/WTS_weight_k_structure/theory'
for k in range(6, 14):
    p = f'{wts}/kern_{k}.pkl'
    print('kern_%d.pkl' % k, time.strftime('%m-%d %H:%M', time.localtime(os.path.getmtime(p))))
with open(wd+'/progress.log','a') as f:
    f.write(f"{datetime.datetime.now():%H:%M:%S} step_13a: verify_W7.py rerun rc={r.returncode} ({(datetime.datetime.now()-t0).seconds}s)\n")
