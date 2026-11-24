
import subprocess, os, time
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
PY = '<project>/.venv/bin/python'
os.chdir(WORK)
with open('step_02.log','w') as log:
    p = subprocess.Popen(['nice','-n','10',PY,'step_02.py'], stdout=log, stderr=subprocess.STDOUT, cwd=WORK)
print("pid", p.pid)
while p.poll() is None:
    time.sleep(25)
    print(open('step_02.log').read()[-1500:])
    print('--- running ---', flush=True)
print("exit", p.returncode)
print(open('step_02.log').read())
