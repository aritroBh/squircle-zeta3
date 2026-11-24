import subprocess, time
PY='<project>/.venv/bin/python'
t0=time.time()
r=subprocess.run(['nice','-n','10',PY,'<project>/research/08_theory/RST1_rest_exponent/verify_rst1.py'],
                 capture_output=True,text=True,timeout=1500,
                 cwd='<project>/research/08_theory/RST1_rest_exponent')
el=time.time()-t0
print(f"rc={r.returncode} elapsed={el:.0f}s")
print(r.stdout)
print("STDERR tail:", r.stderr[-1500:])
open('verify_rst1_rerun.log','w').write(r.stdout+"\nSTDERR:\n"+r.stderr)
