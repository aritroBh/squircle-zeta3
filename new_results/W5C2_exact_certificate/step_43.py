import subprocess, os, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
# robust background launcher: run sectors sequentially under one python driver via shell loop
out = open('step_01.out','wb')
cmd = ['/bin/bash','-c','for s in pp pm mp mm; do <project>/.venv/bin/python step_01.py $s || exit 1; done; echo ALL_SECTORS_DONE']
p = subprocess.Popen(['nice','-n','10'] + cmd, cwd=os.getcwd(), stdout=out, stderr=subprocess.STDOUT,
                     start_new_session=True)
print('launched pid', p.pid)
open('step_01.pid','w').write(str(p.pid))
time.sleep(30)
print(open('step_01.out').read())