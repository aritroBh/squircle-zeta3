import subprocess, datetime, time, os
work = '<project>/research/08_theory/TOW3R_review/'
PY='<project>/.venv/bin/python'
out=''
for args in [('prove_adj34.py','3','0'),('prove_adj34.py','3','1'),('prove_adj34.py','4','2')]:
    t0=time.time()
    r = subprocess.run(['nice','-n','10',PY]+list(args),capture_output=True,text=True,cwd=work+'tow3_copy',timeout=1700)
    line=f'{args}: rc={r.returncode} out={r.stdout.strip()[-200:]} err={r.stderr.strip()[-150:]} ({time.time()-t0:.0f}s)'
    print(line); out+=line+'\n'
open(work+'step_16a.out','w').write(out)
with open(work+'progress.log','a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} step_16a rerun prove_adj34 (3,0),(3,1),(4,2) from TOW3 copies\n")
