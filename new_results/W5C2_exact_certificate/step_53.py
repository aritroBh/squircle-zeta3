import os, time
pid = int(open('<project>/research/08_theory/W5C2_exact_certificate/step_01.pid').read())
try:
    os.kill(pid, 0)
    print('process', pid, 'ALIVE')
except OSError as e:
    print('process', pid, 'DEAD:', e)
# also check the sleeper pid if still around
import subprocess
print(subprocess.run(['bash','-c','tail -2 <project>/research/08_theory/W5C2_exact_certificate/sleeper.log'], capture_output=True, text=True).stdout)