import subprocess, time
PY='<project>/.venv/bin/python'
WORK='<project>/research/08_theory/A3h2_w4_missing'
t=time.time()
for p,tag in [(1000000009,'b'),(1000000033,'b')]:
    r=subprocess.run([PY,WORK+'/rref_driver.py',str(p),tag],capture_output=True,text=True,timeout=1200)
    print(r.stdout, r.stderr[-500:] if r.returncode else '')
print('total',time.time()-t,'s')
