import os, subprocess, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
out = open('verify_W5_final.out','wb')
p = subprocess.Popen(['<project>/.venv/bin/python','verify_W5.py'],
                     cwd=os.getcwd(), stdout=out, stderr=subprocess.STDOUT, start_new_session=True)
print('launched final verify, pid', p.pid)
time.sleep(90)
print(open('verify_W5_final.out').read())
print(subprocess.run(['bash','-c','ls -la *.py *.md *.json 2>/dev/null | head -40'], capture_output=True, text=True).stdout)