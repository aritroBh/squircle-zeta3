import os, shutil, subprocess, datetime, time
base = '<project>/research/08_theory/'
work = base+'TOW3R_review/'
# TOW2R's tow2_copy already has TOW2's files; check
c = base+'TOW2R_review/tow2_copy/'
print(os.path.exists(c), [f for f in ['step_48.py','step_49.py','step_52.py','tow2lib.py'] if os.path.exists(c+f)] if os.path.exists(c) else '')
PY='<project>/.venv/bin/python'
out=''
for args in [('step_48.py',), ('step_52.py',)]:
    t0=time.time()
    r = subprocess.run(['nice','-n','10',PY]+list(args),capture_output=True,text=True,cwd=c,timeout=1700)
    line=f'TOW2R copy {args}: rc={r.returncode} out={r.stdout.strip()[-200:]} err={r.stderr.strip()[-150:]} ({time.time()-t0:.0f}s)'
    print(line); out+=line+'\n'
open(work+'step_17.out','w').write(out)
with open(work+'progress.log','a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} step_17 rerun TOW2 V(3) and slope(4,1) from TOW2R tow2_copy\n")
