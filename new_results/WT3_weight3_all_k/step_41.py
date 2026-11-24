
import subprocess, os, time
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
PY = '<project>/.venv/bin/python'
os.chdir(WORK)
with open('step_04.log','w') as log:
    p = subprocess.Popen(['nice','-n','10',PY,'step_04.py'], stdout=log, stderr=subprocess.STDOUT, cwd=WORK)
print("pid", p.pid)
while p.poll() is None:
    time.sleep(30)
    print(open('step_04.log').read()[-600:]); print('--- running ---', flush=True)
print("exit", p.returncode)
out = open('step_04.log').read()
print(out[-4000:])
