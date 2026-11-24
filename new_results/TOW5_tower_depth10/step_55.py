import subprocess, time, os
work = '<project>/research/08_theory/TOW5_tower_depth10'
PY = '<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run(['nice','-n','10',PY,'forced910.py','10'], cwd=work, capture_output=True, text=True, timeout=1700)
print(r.stdout)
print("ERR", r.stderr[-800:] if r.stderr else "")
print("elapsed", round(time.time()-t0,1), "rc", r.returncode)
with open(work+'/progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: forced910.py 10 rc={r.returncode} ({round(time.time()-t0,1)}s)\n")
