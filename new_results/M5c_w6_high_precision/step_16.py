print(open('m5b_indep_check.py').read()[:6000])
import subprocess
print(subprocess.run(['uptime'],capture_output=True,text=True).stdout)
