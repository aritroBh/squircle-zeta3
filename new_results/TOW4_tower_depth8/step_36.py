import subprocess, os, time
work = '<project>/research/08_theory/TOW4_tower_depth8'
PY = '<project>/.venv/bin/python'
# remove stale pkl so we start fresh
if os.path.exists(work + '/forced_F78.pkl'): os.remove(work + '/forced_F78.pkl')
logf = open(work + '/forced78.log','w')
proc = subprocess.Popen(['nice','-n','10', PY, work + '/forced78.py'], stdout=logf, stderr=subprocess.STDOUT, cwd=work)
open(work + '/forced78.pid','w').write(str(proc.pid))
print("launched pid", proc.pid)
time.sleep(20)
print(open(work + '/forced78.log').read())
