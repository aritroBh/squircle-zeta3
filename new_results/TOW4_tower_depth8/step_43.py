import subprocess, os, time
work = '<project>/research/08_theory/TOW4_tower_depth8'
PY = '<project>/.venv/bin/python'
if os.path.exists(work + '/forced_F78.pkl'): os.remove(work + '/forced_F78.pkl')
logf = open(work + '/forced78.log','w')
proc = subprocess.Popen(['nice','-n','10', PY, work + '/forced78.py'], stdout=logf, stderr=subprocess.STDOUT, cwd=work)
open(work + '/forced78.pid','w').write(str(proc.pid))
print("launched pid", proc.pid)
# poll up to 200s
for _ in range(40):
    time.sleep(5)
    txt = open(work + '/forced78.log').read()
    if 'VERDICT' in txt or 'Traceback' in txt: break
print(txt[txt.find('delta=7 done'):])
