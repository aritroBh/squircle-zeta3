import os, signal, subprocess, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
try:
    os.kill(int(open('step_01.pid').read()), signal.SIGKILL); print('killed loop')
except Exception as e: print(e)
time.sleep(2)
out = open('step_01b.out','wb')
p = subprocess.Popen(['nice','-n','10','/bin/bash','-c','for s in pp pm mp mm; do <project>/.venv/bin/python step_01b.py $s || exit 1; done; echo ALL_SECTORS_DONE'],
                     cwd=os.getcwd(), stdout=out, stderr=subprocess.STDOUT, start_new_session=True)
open('step_01b.pid','w').write(str(p.pid))
print('launched step_01b all sectors, pid', p.pid)
time.sleep(30)
print(open('step_01b.out').read())