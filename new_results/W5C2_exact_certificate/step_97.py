import os, subprocess, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
with open('progress.log','a') as f:
    f.write('2026-10-29 04:38:46 step_02 ALL DONE (4 sectors, all checks passed mod p2).\n')
out = open('step_03_test.out','wb')
p = subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','step_03_test.py'],
                     cwd=os.getcwd(), stdout=out, stderr=subprocess.STDOUT, start_new_session=True)
open('step_03_test.pid','w').write(str(p.pid))
print('launched test, pid', p.pid)
time.sleep(200)
print(open('step_03_test.out').read())