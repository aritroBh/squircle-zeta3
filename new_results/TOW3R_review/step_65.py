import subprocess, datetime, time
base = '<project>/research/08_theory/'
work = base+'TOW3R_review/'
PY='<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run(['nice','-n','10',PY,'step_49.py'],capture_output=True,text=True,cwd=base+'TOW2R_review/tow2_copy',timeout=1700)
print(f'V(4): rc={r.returncode} out={r.stdout.strip()[-200:]} err={r.stderr.strip()[-150:]} ({time.time()-t0:.0f}s)')
open(work+'step_18.out','w').write(f'V(4) rerun: rc={r.returncode} {r.stdout.strip()[-200:]} ({time.time()-t0:.0f}s)')
with open(work+'progress.log','a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} step_18 rerun TOW2 V(4) step_49\n")
