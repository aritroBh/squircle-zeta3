
import subprocess, os, time
Z3DIR = 'Zeta' + chr(40) + '3' + chr(41)
BASE = os.path.join('<home>/Downloads/Research', Z3DIR)
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
PY = os.path.join(BASE, '.venv/bin/python')
os.chdir(WORK)
with open('verify_wt3.log','w') as log:
    p = subprocess.Popen(['nice','-n','10',PY,'verify_wt3.py'], stdout=log, stderr=subprocess.STDOUT, cwd=WORK)
print("pid", p.pid)
while p.poll() is None:
    time.sleep(30)
    print(open('verify_wt3.log').read()[-500:]); print('--- running ---', flush=True)
print("exit", p.returncode)
print(open('verify_wt3.log').read())
