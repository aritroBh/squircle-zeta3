import subprocess, os, time
PY = '<project>/.venv/bin/python'
logf = open('resultant.log','w')
p = subprocess.Popen(['nice','-n','10',PY,'run_resultant.py'], stdout=logf, stderr=subprocess.STDOUT,
                     cwd='<project>/research/08_theory/CHN3_p5_locus')
print("PID:", p.pid)
time.sleep(20)
print(open('resultant.log').read())
