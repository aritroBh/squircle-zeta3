import subprocess, datetime, time
work = '<project>/research/08_theory/TOW3R_review/'
PY='<project>/.venv/bin/python'
out_all=''
for (delta, s) in [(5,0),(5,1),(5,2)]:
    t0=time.time()
    r = subprocess.run(['nice','-n','10',PY,'prove_adj56.py',str(delta),str(s)],
                       capture_output=True,text=True,cwd=work+'tow3_copy',timeout=1700)
    line = f'adj({delta},{s}): rc={r.returncode} out={r.stdout.strip()} err={r.stderr.strip()[-200:]} ({time.time()-t0:.0f}s)'
    print(line); out_all+=line+'\n'
open(work+'step_09a.out','w').write(out_all)
with open(work+'progress.log','a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} step_09a rerun prove_adj56 delta=5\n")
