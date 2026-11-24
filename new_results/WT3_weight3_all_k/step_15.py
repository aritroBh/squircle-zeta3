
import subprocess, os, time
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
PY = '<project>/.venv/bin/python'
os.chdir(WORK)
# launch step_01 under nice in background
with open('step_01.log','w') as log:
    p = subprocess.Popen(['nice','-n','10',PY,'step_01.py'], stdout=log, stderr=subprocess.STDOUT, cwd=WORK)
print("pid", p.pid)
# poll
while p.poll() is None:
    time.sleep(20)
    print(open('step_01.log').read()[-800:])
    print('--- still running ---', flush=True)
print("exit", p.returncode)
print(open('step_01.log').read())
