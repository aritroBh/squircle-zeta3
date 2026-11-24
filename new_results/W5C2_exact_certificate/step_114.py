import os, subprocess, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
with open('progress.log','a') as f:
    f.write('2026-10-29 06:05:53 step_03 ALL SOLVED (4 sectors, exact Dixon, heights ~1990-2081 digits, all cross-checks mod p2 OK). Total solve ~54 min.\n')
out = open('step_04a.out','wb')
p = subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','step_04a.py'],
                     cwd=os.getcwd(), stdout=out, stderr=subprocess.STDOUT, start_new_session=True)
open('step_04a.pid','w').write(str(p.pid))
print('launched step_04a, pid', p.pid)
time.sleep(150)
print(open('step_04a.out').read())