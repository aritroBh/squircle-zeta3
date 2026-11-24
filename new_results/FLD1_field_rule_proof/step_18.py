import subprocess, time
VENV='<project>/.venv/bin/python'
G1='<project>/research/10_frontier/G1_zeta3_any_angle/'
# Test flat_exact timing for k=5, N=12 first (fast)
t=time.time()
r = subprocess.run(['nice',VENV,'flat_exact.py','12','5'],capture_output=True,text=True,cwd=G1)
print("k<=5,N=12 time", round(time.time()-t,1),"s")
print(r.stdout[-1500:])
print(r.stderr[-500:] if r.returncode else '')