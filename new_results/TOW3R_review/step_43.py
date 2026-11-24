import subprocess, datetime, time
work = '<project>/research/08_theory/TOW3R_review/'
PY='<project>/.venv/bin/python'
out_all=''
for delta in [5,6]:
    t0=time.time()
    r = subprocess.run(['nice','-n','10',PY,'prove_E3b56.py',str(delta)],
                       capture_output=True,text=True,cwd=work+'tow3_copy',timeout=1700)
    line = f'E3({delta}): rc={r.returncode} out={r.stdout.strip()} err={r.stderr.strip()[-200:]} ({time.time()-t0:.0f}s)'
    print(line); out_all+=line+'\n'
open(work+'step_10.out','w').write(out_all)
with open(work+'progress.log','a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} step_10 rerun prove_E3b56 delta=5,6\n")
