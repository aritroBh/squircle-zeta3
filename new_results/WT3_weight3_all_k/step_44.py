
import subprocess, os, time
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
PY = '<project>/.venv/bin/python'
os.chdir(WORK)
with open('step_05.log','w') as log:
    p = subprocess.Popen(['nice','-n','10',PY,'step_05.py'], stdout=log, stderr=subprocess.STDOUT, cwd=WORK)
while p.poll() is None:
    time.sleep(20)
print("exit", p.returncode)
print(open('step_05.log').read())
