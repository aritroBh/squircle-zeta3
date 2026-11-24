import subprocess, time, os
work = '<project>/research/08_theory/TOW5_tower_depth10'
PY = '<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run(['nice','-n','10',PY,'verify_tow5.py'], cwd=work, capture_output=True, text=True, timeout=2400)
print(r.stdout)
print("STDERR tail:", r.stderr[-1500:] if r.stderr else "")
print("elapsed", round(time.time()-t0,1), "rc", r.returncode)
# save log
open(work+'/verify_tow5.log','w').write(r.stdout + "\nSTDERR:\n" + (r.stderr or ""))
with open(work+'/progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: verify_tow5.py rc={r.returncode} ({round(time.time()-t0,1)}s)\n")
