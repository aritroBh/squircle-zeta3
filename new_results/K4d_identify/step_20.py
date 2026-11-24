import subprocess, os, time, sys
PY = '<home>)/.venv/bin/python'
WORK = 'research/08_theory/K4d_identify'
os.chdir(WORK)

# Test detached subprocess survival: launch a 6s sleeper that writes a file
sleeper = os.path.join(WORK,'_sleeper.py')
open(sleeper,'w').write(
"import time,os\n"
"f=open('%s','w'); f.write('start'); f.close()\n"%os.path.join(WORK,'_sleep_out.txt')+
"time.sleep(6)\n"
"f=open('%s','w'); f.write('done'); f.close()\n"%os.path.join(WORK,'_sleep_out.txt'))
if os.path.exists('_sleep_out.txt'): os.remove('_sleep_out.txt')
p = subprocess.Popen([PY,'_sleeper.py'], cwd=WORK,
    stdout=open('_sleep_log.txt','w'), stderr=subprocess.STDOUT,
    start_new_session=True)
print("launched pid", p.pid, "returning immediately")
print("out exists now:", os.path.exists('_sleep_out.txt'))
