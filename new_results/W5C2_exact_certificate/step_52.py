import subprocess, os, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
st = os.stat('step_01.out'); print('out mtime:', time.ctime(st.st_mtime))
r = subprocess.run(['bash','-c','top -l 1 -n 5 2>&1 | head -20'], capture_output=True, text=True)
print(r.stdout[:800])