import subprocess, os, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
# clean any stale state; test sleeper survival across turn boundary
open('sleeper.log','w').write('start\n')
code = "import time\nfor i in range(12):\n    time.sleep(10)\n    open('sleeper.log','a').write(f'tick {i} {time.ctime()}\\n')\n"
open('sleeper.py','w').write(code)
p = subprocess.Popen(['<project>/.venv/bin/python','sleeper.py'],
                     cwd=os.getcwd(), stdout=open('sleeper.out','wb'), stderr=subprocess.STDOUT,
                     start_new_session=True)
print('launched pid', p.pid)