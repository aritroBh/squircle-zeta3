import subprocess, datetime, time, os
work = '<project>/research/08_theory/TOW3R_review/'
base = '<project>/research/08_theory/'
print(open(work+'tow3_copy/prove_adj34.py').read()[:400])
PY='<project>/.venv/bin/python'
out=''
for args in [('prove_adj34.py','4','0'),('prove_adj34.py','4','1'),('step_48.py',),('step_52.py',)]:
    t0=time.time()
    r = subprocess.run(['nice','-n','10',PY]+list(args),capture_output=True,text=True,cwd=work+'tow3_copy',timeout=1700)
    line=f'{args}: rc={r.returncode} out={r.stdout.strip()[-250:]} err={r.stderr.strip()[-150:]} ({time.time()-t0:.0f}s)'
    print(line); out+=line+'\n'
open(work+'step_16b.out','w').write(out)
with open(work+'progress.log','a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} step_16b rerun adj(4,0),(4,1), V3(step_48), slope41(step_52)\n")
