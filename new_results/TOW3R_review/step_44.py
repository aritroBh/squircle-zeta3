import shutil, subprocess, datetime, time
base = '<project>/research/08_theory/'
work = base+'TOW3R_review/'
for f in ['prove_S12.py','prove_S34.py','prove_tops.py']:
    shutil.copy(base+'TOW3_tower_depth6/'+f, work+'tow3_copy/'+f)
PY='<project>/.venv/bin/python'
out_all=''
for args in [('prove_S12.py',),('prove_S34.py','3'),('prove_S34.py','4'),('prove_tops.py',)]:
    t0=time.time()
    r = subprocess.run(['nice','-n','10',PY]+list(args), capture_output=True,text=True,cwd=work+'tow3_copy',timeout=1700)
    line = f'{args}: rc={r.returncode} out={r.stdout.strip()[-300:]} err={r.stderr.strip()[-200:]} ({time.time()-t0:.0f}s)'
    print(line); out_all+=line+'\n'
open(work+'step_11.out','w').write(out_all)
with open(work+'progress.log','a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} step_11 rerun TOW2 base certs S12/S34/tops from TOW3 copies\n")
