import subprocess, os, time
work = '<project>/research/08_theory/TOW4_tower_depth8'
PY = '<project>/.venv/bin/python'
logf = open(work + '/verify_tow4.log','w')
proc = subprocess.Popen(['nice','-n','10', PY, work + '/verify_tow4.py'], stdout=logf, stderr=subprocess.STDOUT, cwd=work)
print("launched", proc.pid)
time.sleep(120)
print(open(work + '/verify_tow4.log').read())
