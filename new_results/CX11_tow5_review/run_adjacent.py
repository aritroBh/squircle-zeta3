"""Run all eleven exact fast certificates, preserving per-case return status."""
import subprocess,sys
from pathlib import Path
with open('report.md','a') as report:
 report.write('\n## Item 4 live execution record\n\n')
for dl in (9,10):
 for s in range(dl//2+1):
  name=f'adj_{dl}_{s}.log'
  with open(name,'w') as out:
   r=subprocess.run([sys.executable,'-B','run_certificate.py','prove_adj910_fast.py',str(dl),str(s)],stdout=out,stderr=subprocess.STDOUT,timeout=300)
  print(dl,s,'exit',r.returncode,flush=True)
  with open('report.md','a') as report: report.write(f'- Fast adjacent ({dl},{s}): exit {r.returncode}; evidence `{name}`.\n')
  assert r.returncode==0,name
print('All 11 adjacent cases PASS.',flush=True)
