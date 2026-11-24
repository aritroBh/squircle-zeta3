import os, subprocess, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
out = open('verify_W5.out','wb')
p = subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','verify_W5.py'],
                     cwd=os.getcwd(), stdout=out, stderr=subprocess.STDOUT, start_new_session=True)
open('verify.pid','w').write(str(p.pid))
print('launched verify_W5, pid', p.pid)
time.sleep(180)
print(open('verify_W5.out').read())