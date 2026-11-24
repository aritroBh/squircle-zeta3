import subprocess, datetime, time
work = '<project>/research/08_theory/TOW3R_review/'
PY='<project>/.venv/bin/python'
out_all = ''
for (delta, s) in [(5,0),(5,1),(5,2),(6,0),(6,1),(6,2),(6,3)]:
    t0=time.time()
    r = subprocess.run(['nice','-n','10',PY,'prove_V56.py',str(delta),str(s)],
                       capture_output=True,text=True,cwd=work+'tow3_copy',timeout=1700)
    line = f'V({delta},{s}): rc={r.returncode} out={r.stdout.strip()} err={r.stderr.strip()[-200:]} ({time.time()-t0:.0f}s)'
    print(line); out_all += line+'\n'
open(work+'step_08.out','w').write(out_all)
with open(work+'progress.log','a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} step_08 rerun prove_V56.py all 7 identities\n")
