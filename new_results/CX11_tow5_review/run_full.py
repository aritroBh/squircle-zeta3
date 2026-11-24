import subprocess,sys
for script in ('check_full910.py','check_full910_ext.py'):
 name=script.removesuffix('.py')+'.log'
 with open(name,'w') as out:
  r=subprocess.run([sys.executable,'-B','run_certificate.py',script],stdout=out,stderr=subprocess.STDOUT,timeout=300)
 print(script,'exit',r.returncode,flush=True)
 assert r.returncode==0
 with open('report.md','a') as f: f.write(f'\nFull-law certificate `{script}` completed with exit 0; see `{name}` for counts.\n')
